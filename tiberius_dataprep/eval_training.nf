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
params.config_yaml = '../config/config.yaml'
params.use_test    = false

def cfg        = new groovy.yaml.YamlSlurper().parse(file(params.config_yaml))
def GENOME_DIR = cfg.genome_dir as String
def WORK_DIR   = cfg.work_dir   as String
def splitName  = params.use_test ? 'test' : 'val'
def speciesList = cfg.species_split?."${splitName}" ?: []

if( !speciesList )
    error "species_split.${splitName} is empty or missing in ${params.config_yaml}"

if( !cfg.training )
    error "training: section is empty or missing in ${params.config_yaml}"

// Derive per-species genome + reference annotation paths from the dataprep
// conventions: genome under genome_dir, ref annot is dataprep's output GTF.
def speciesMeta = speciesList.collect { sp ->
    [
        sp,
        "${GENOME_DIR}/${sp}.genome.fa".toString(),
        "${WORK_DIR}/annot_gtf/${sp}_longest.gtf".toString()
    ]
}

// epochInfoList rows: [ trainName, idx, epochDir (str), evalDir (str), hmmFlag ]
def epochInfoList = []
cfg.training.each { trainName, trainAttrs ->

    def hmmFlag = (trainAttrs.hmm ?: false)
    def evalDir = trainAttrs.eval_dir?.toString()
    if( !evalDir )
        error "training.${trainName}.eval_dir is missing in ${params.config_yaml}"

    if( params.use_test ) {
        (trainAttrs.test_weights ?: []).eachWithIndex { twPath, i ->
            epochInfoList << [
                trainName,
                "test_${i}",
                file(twPath).toString(),
                evalDir,
                hmmFlag
            ]
        }
        return  // skip the weights_dirs scan for this trainName
    }

    (trainAttrs.weights_dirs ?: []).eachWithIndex { weightsDirPath, idx ->
        def baseDir = file(weightsDirPath.toString())
        if( baseDir.isDirectory() ) {
            baseDir.listFiles()
                .findAll { it.isDirectory() && it.name.startsWith('epoch_') }
                .each { epochSubdir ->
                    epochInfoList << [
                        trainName,
                        idx,
                        epochSubdir.toString(),
                        evalDir,
                        hmmFlag
                    ]
                }
        }
    }
}

// Channel of (species, genome.fa, annot.gtf) for the chosen split.
Channel
    .from( speciesMeta )
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
