#!/usr/bin/env nextflow
/*
 * Requires: nextflow >= 23.04
 */

nextflow.enable.dsl = 2

params.container   = 'docker://larsgabriel23/tiberius:dev'
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
def epochInfoList = []
yamlConfig.Training.each { trainName, trainAttrs ->

    def hmmFlag = (trainAttrs.HMM ?: false)

    if( params.use_test && trainAttrs.TestWeights ) {
        trainAttrs.TestWeights.each { twKey, twPathObj ->

            def epochSubdir  = file(twPathObj).toString()               // epoch_* dir
            def weightsDir = file(twPathObj).parent.toString()  // its parent
            
            epochInfoList << [
                trainName,
                twKey,              // we’ll use the YAML key as “idx”
                epochSubdir,             // epochDir
                trainAttrs.EvalDir.toString(),
                weightsDir,         // weightsDirPath
                hmmFlag
            ]
        }
        return                                   // go to next trainName
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
                        weightsDirPath.toString(), 
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

// Channel of (trainName, idx, epochDirAsFile, evalDir)
Channel
    .from(epochInfoList)
    .map { trainEpochTuple ->
        def (trainName, idx, epochPath, evalDir, weightsDirPath, hmmFlag) = trainEpochTuple
        tuple(trainName, idx, file(epochPath), file(evalDir), file(weightsDirPath), hmmFlag)
    }
    .set { EPOCH_INFO }

// Cross-join species metadata with every (trainName, epochDir)
SPECIES_META
    .combine( EPOCH_INFO )                        
    .map { speciesName, genomeFa, annotGtf,        
           trainName,  idx, epochDir, evalDir,
           weightsDirPath, hmmFlag ->

        tuple(
            trainName, idx,
            epochDir, evalDir, weightsDirPath, hmmFlag,
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

    publishDir {"${evalDir}/"}, mode:'copy', pattern: '*.gtf'

    container params.container
    
    memory '190 GB'
    input:
        tuple(
            val(trainName),  
            val(idx), 
            path(epochDir), 
            val(evalDir),
            path(weightsDirPath),
            val(hmmFlag),
            val(speciesName),
            path(genomeFa),
            path(annotGtf)
        )

    output:
        // emit the metadata *together with* the resulting GTF
        tuple (
            path("${idx}.${epochDir}.${speciesName}.gtf"),
            val(trainName),  
            val(idx), 
            path(epochDir), 
            val(evalDir),
            path(weightsDirPath),
            val(speciesName),
            path(annotGtf)
        )

    script:
    """
    tiberius.py \
        --genome ${genomeFa} \
        ${hmmFlag ? '--model_old' : '--model_lstm_old'} ${weightsDirPath}/${epochDir} \
        --out ${idx}.${epochDir}.${speciesName}.gtf
    """
}

//--------------------------------------------------------------------------
// PROCESS 2 ─ RUN_GFFCOMPARE
//--------------------------------------------------------------------------
process RUN_GFFCOMPARE {

    publishDir {"${evalDir}/"}, mode:'copy', pattern: '*.stats'

    input:
        tuple (
            path(gtf),
            val(trainName),  
            val(idx), 
            path(epochDir), 
            val(evalDir),
            path(weightsDirPath),
            val(speciesName),
            path(annotGtf)
        )

    output:
        path "${idx}.${epochDir}.${speciesName}.stats"

    script:
    """
    grep CDS ${annotGtf} > annot_cds.gtf
    gffcompare \
      -r annot_cds.gtf \
      --strict-match -e 0 -T --no-merge \
      ${gtf} \
      -o ${idx}.${epochDir}.${speciesName}.stats
    """
}

workflow {
    JOB_INPUTS | RUN_TIBERIUS | RUN_GFFCOMPARE
}