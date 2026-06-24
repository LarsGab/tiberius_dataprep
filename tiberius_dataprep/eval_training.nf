#!/usr/bin/env nextflow
/*
 * Requires: nextflow >= 23.04
 */

nextflow.enable.dsl = 2

// Supported Tiberius versions -> pinned container tag.
// Override the whole image with --container docker://my/image:tag if needed.
def TIBERIUS_TAG_FOR = [ '1.x': '1.1.8', '2.x': '2.0.6' ]
params.tiberius_version = '2.x'
if (!TIBERIUS_TAG_FOR.containsKey(params.tiberius_version)) {
    error "params.tiberius_version must be one of ${TIBERIUS_TAG_FOR.keySet()} (got '${params.tiberius_version}')"
}
params.container   = "docker://larsgabriel23/tiberius:${TIBERIUS_TAG_FOR[params.tiberius_version]}"
params.config      = '../config/config_train_eval.yaml'

import groovy.yaml.YamlSlurper

def yamlConfig = new YamlSlurper().parse(file(params.config))
def dataSpecies = params.use_test ? yamlConfig.Test_Data
                                  : yamlConfig.Val_Data

if( !dataSpecies )
    throw new IllegalArgumentException(
        "Requested ${params.use_test ? 'Test_Data' : 'Val_Data'} section " +
        "is missing in config_train_eval.yaml" )

params.species_meta = dataSpecies.collect { speciesName, attrs ->
    [
        speciesName,
        attrs.Genome.toString(),
        attrs.RefAnnot.toString()
    ]
}

// 3) Scan all Training → WeightsDir → epoch_* subdirs to collect every epoch directory
// epochInfoList rows: [ trainName, idx, epochDir (str), evalDir (str), hmmFlag ]
def epochInfoList = []
yamlConfig.Training.each { trainName, trainAttrs ->

    def hmmFlag = (trainAttrs.HMM ?: false)

    if( params.use_test && trainAttrs.TestWeights ) {
        trainAttrs.TestWeights.each { twKey, twPathObj ->
            epochInfoList << [
                trainName,
                twKey,                              // YAML key used as idx
                file(twPathObj).toString(),
                trainAttrs.EvalDir.toString(),
                hmmFlag
            ]
        }
        return  // skip the WeightsDir loop for this trainName
    }

    trainAttrs.WeightsDir.each { idx, weightsDirPath ->
        def baseDir = file(weightsDirPath.toString())
        if( baseDir.isDirectory() ) {
            baseDir.listFiles()
                .findAll { it.isDirectory() && it.name.startsWith('epoch_') }
                .each { epochSubdir ->
                    epochInfoList << [
                        trainName,
                        idx,
                        epochSubdir.toString(),
                        trainAttrs.EvalDir.toString(),
                        hmmFlag
                    ]
                }
        }
    }
}

// Channel of (species, genome.fa, annot.gtf) using Val_Data only
Channel
    .from( params.species_meta )
    .map { row -> tuple( row[0], file(row[1]), file(row[2]) ) }
    .set { SPECIES_META }

// Channel of (trainName, idx, epochDir, evalDir, hmmFlag)
Channel
    .from(epochInfoList)
    .map { trainName, idx, epochPath, evalDir, hmmFlag ->
        tuple(trainName, idx, file(epochPath), evalDir, hmmFlag)
    }
    .set { EPOCH_INFO }

// Cross-join species metadata with every (trainName, epochDir)
SPECIES_META
    .combine( EPOCH_INFO )
    .map { speciesName, genomeFa, annotGtf,
           trainName, idx, epochDir, evalDir, hmmFlag ->

        tuple(
            trainName, idx,
            epochDir, evalDir, hmmFlag,
            speciesName, genomeFa, annotGtf
        )
    }
    .set { JOB_INPUTS }

//--------------------------------------------------------------------------
// Process 1: run_tiberius  ≈  rule run_tiberius
//--------------------------------------------------------------------------
// ─────────── 4-A.  Process RUN_TIBERIUS – add species & genome.fa ───────
process RUN_TIBERIUS {
    label 'gpu'

    publishDir { "${evalDir}/" }, mode: 'copy', pattern: '*.gtf'

    container params.container
    storeDir { "cache/${task.process}/${trainName}/${speciesName}/${idx}/${epochDir.name}/" }
    memory '190 GB'
    input:
        tuple(
            val(trainName),
            val(idx),
            path(epochDir),
            val(evalDir),
            val(hmmFlag),
            val(speciesName),
            path(genomeFa),
            path(annotGtf)
        )

    output:
        tuple(
            path("${idx}.${epochDir.name}.${speciesName}.gtf"),
            val(trainName),
            val(idx),
            path(epochDir),
            val(evalDir),
            val(speciesName),
            path(annotGtf)
        )

    script:
    """
    tiberius.py \\
        --genome ${genomeFa} \\
        ${hmmFlag ? '--model_old' : '--model_lstm_old'} ${epochDir} \\
        --out ${idx}.${epochDir.name}.${speciesName}.gtf
    """
}

//--------------------------------------------------------------------------
// PROCESS 2 ─ RUN_GFFCOMPARE
//--------------------------------------------------------------------------
process RUN_GFFCOMPARE {

    publishDir { "${evalDir}/" }, mode: 'copy', pattern: '*.stats'
    storeDir { "cache/${task.process}/${trainName}/${speciesName}/${idx}/${epochDir.name}/" }

    input:
        tuple(
            path(gtf),
            val(trainName),
            val(idx),
            path(epochDir),
            val(evalDir),
            val(speciesName),
            path(annotGtf)
        )

    output:
        path "${idx}.${epochDir.name}.${speciesName}.stats"

    script:
    """
    grep CDS ${annotGtf} > annot_cds.gtf
    gffcompare \\
      -r annot_cds.gtf \\
      --strict-match -e 0 -T --no-merge \\
      ${gtf} \\
      -o ${idx}.${epochDir.name}.${speciesName}.stats
    """
}

workflow {
    JOB_INPUTS | RUN_TIBERIUS | RUN_GFFCOMPARE
}