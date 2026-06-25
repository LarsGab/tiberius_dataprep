#!/usr/bin/env nextflow
/*
 * Requires: nextflow >= 23.04
 */

nextflow.enable.dsl = 2

////////////////////////////////////////////////////////////////////////
//                     PARAMETER DEFAULTS                              //
////////////////////////////////////////////////////////////////////////

// Override the whole image with --container docker://my/image:tag if needed.
params.tiberius_version = '2.x'
// Accept the legacy --config flag (eval_training only) and --configYAML
// (dataprep) as aliases for --config_yaml. "First assignment wins" in
// Nextflow 26, so a CLI --config_yaml beats the alias.
params.config_yaml      = params.containsKey('config')     ? params.config
                        : params.containsKey('configYAML') ? params.configYAML
                        : '../config/config.yaml'
params.use_test         = false
// NB: do NOT declare params.container at top level. In Nextflow 26 the first
// assignment wins, so a null default would shadow the workflow's resolution.

////////////////////////////////////////////////////////////////////////
//                            HELPERS                                  //
////////////////////////////////////////////////////////////////////////

def tiberiusImage(String version) {
    def tagFor = [ '1.x': '1.1.8', '2.x': '2.0.6' ]
    if (!tagFor.containsKey(version))
        error "params.tiberius_version must be one of ${tagFor.keySet()} (got '${version}')"
    return "docker://larsgabriel23/tiberius:${tagFor[version]}"
}

////////////////////////////////////////////////////////////////////////
//                            PROCESSES                               //
////////////////////////////////////////////////////////////////////////

process RUN_TIBERIUS {
    label 'gpu'

    publishDir { "${evalDir}/" }, mode: 'copy', pattern: '*.gtf'

    container { params.container }
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

process RUN_GFFCOMPARE {

    publishDir { "${evalDir}/" }, mode: 'copy', pattern: '*.stats'

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

////////////////////////////////////////////////////////////////////////
//                            WORKFLOW                                //
////////////////////////////////////////////////////////////////////////

workflow {

    // Resolve container from tiberius_version unless the user pinned it.
    if (!params.container) {
        params.container = tiberiusImage(params.tiberius_version)
    }

    def cfg        = new groovy.yaml.YamlSlurper().parse(file(params.config_yaml))
    def GENOME_DIR = cfg.genome_dir as String
    def WORK_DIR   = cfg.work_dir   as String
    def splitName  = params.use_test ? 'test' : 'val'
    def speciesList = cfg.species_split?."${splitName}" ?: []

    if( !speciesList )
        error "species_split.${splitName} is empty or missing in ${params.config_yaml}"

    if( !cfg.training )
        error "training: section is empty or missing in ${params.config_yaml}"

    // Derive per-species genome + reference annotation from dataprep conventions.
    def speciesMeta = speciesList.collect { sp ->
        def genomeCandidates = [
            file("${GENOME_DIR}/${sp}.genome.fa"),
            file("${GENOME_DIR}/${sp}.fa"),
            file("${GENOME_DIR}/${sp}.fasta"),
        ]
        def genomeFile = genomeCandidates.find { it.exists() }
        if (!genomeFile)
            error "No genome FASTA (.genome.fa, .fa or .fasta) found for ${sp} in ${GENOME_DIR}"

        [
            sp,
            genomeFile.toString(),
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

    def speciesMetaCh = Channel
        .from(speciesMeta)
        .map { row -> tuple( row[0], file(row[1]), file(row[2]) ) }

    def epochInfoCh = Channel
        .from(epochInfoList)
        .map { trainName, idx, epochPath, evalDir, hmmFlag ->
            tuple(trainName, idx, file(epochPath), evalDir, hmmFlag)
        }

    def jobInputs = speciesMetaCh
        .combine( epochInfoCh )
        .map { speciesName, genomeFa, annotGtf,
               trainName, idx, epochDir, evalDir, hmmFlag ->

            tuple(
                trainName, idx,
                epochDir, evalDir, hmmFlag,
                speciesName, genomeFa, annotGtf
            )
        }

    jobInputs | RUN_TIBERIUS | RUN_GFFCOMPARE
}
