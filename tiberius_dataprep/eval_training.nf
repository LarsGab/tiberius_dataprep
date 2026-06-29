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

// Probe an epoch directory and return:
//   'new' -- new layout with model_config.json + weights.h5 / model.weights.h5
//   'old' -- legacy SavedModel layout with saved_model.pb + variables/
// Errors out if neither marker is present.
def detectModelFormat(epochDir) {
    def dir = epochDir instanceof File ? epochDir : new File(epochDir.toString())
    if (new File(dir, 'model_config.json').exists()) return 'new'
    if (new File(dir, 'saved_model.pb').exists())    return 'old'
    error "Cannot determine model format for ${dir}: expected either 'model_config.json' (new layout) or 'saved_model.pb' (old SavedModel layout)."
}

////////////////////////////////////////////////////////////////////////
//                            PROCESSES                               //
////////////////////////////////////////////////////////////////////////

process RUN_TIBERIUS {
    label 'gpu'

    // The pattern matches our 4-segment output filename (idx.epoch.species.gtf)
    // but not the 2-segment staged reference annotation (<species>_longest.gtf).
    publishDir { "${evalDir}/" }, mode: 'copy', pattern: '*.*.*.gtf'

    container { params.container }
    memory '190 GB'
    input:
        tuple(
            val(trainName),
            val(idx),
            path(epochDir),
            val(evalDir),
            val(hmmFlag),
            val(modelFormat),
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
    // Pick the Tiberius flag from the epoch directory layout (see
    // detectModelFormat). HMM presence in 'new' layouts is encoded in
    // model_config.json, so the hmm flag only matters for legacy dirs.
    def modelFlag
    if (modelFormat == 'new') {
        modelFlag = '--model'
    } else {
        modelFlag = hmmFlag ? '--model_old' : '--model_lstm_old'
    }
    """
    tiberius.py \\
        --genome ${genomeFa} \\
        ${modelFlag} ${epochDir} \\
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
        tuple(
            val(trainName),
            val(evalDir),
            path("${idx}.${epochDir.name}.${speciesName}.stats")
        )

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

process PLOT_EVAL {
    tag { "${trainName}_${mode}" }
    // No container directive: matplotlib + pyyaml are not in the Tiberius
    // image, so this process runs on the host Python env. Install with
    // `pip install -e .` from the repo root before running the pipeline.
    publishDir { "${evalDir}/" }, mode: 'copy', pattern: '*.png'

    cpus   1
    memory '2 GB'

    input:
        tuple(
            val(trainName),
            val(evalDir),
            val(speciesCsv),
            val(mode),
            path(stats)
        )

    output:
        path "*.png"

    script:
    """
    plot_eval.py \\
        --stats_dir . \\
        --species "${speciesCsv}" \\
        --mode ${mode} \\
        --title "${trainName}" \\
        --out ${trainName}_${mode}_acc.png
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
    def speciesList = (cfg.species_split ?: [:])[splitName] ?: []

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

    // epochInfoList rows: [ trainName, idx, epochDir, evalDir, hmmFlag, modelFormat ]
    // modelFormat is 'new' or 'old', auto-detected from the directory contents.
    def epochInfoList = []
    cfg.training.each { trainName, trainAttrs ->

        def hmmFlag = (trainAttrs.hmm ?: false)
        def evalDir = trainAttrs.eval_dir?.toString()
        if( !evalDir )
            error "training.${trainName}.eval_dir is missing in ${params.config_yaml}"

        if( params.use_test ) {
            (trainAttrs.test_weights ?: []).eachWithIndex { twPath, i ->
                def epochDir = file(twPath)
                epochInfoList << [
                    trainName,
                    "test_${i}",
                    epochDir.toString(),
                    evalDir,
                    hmmFlag,
                    detectModelFormat(epochDir)
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
                            hmmFlag,
                            detectModelFormat(epochSubdir)
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
        .map { trainName, idx, epochPath, evalDir, hmmFlag, modelFormat ->
            tuple(trainName, idx, file(epochPath), evalDir, hmmFlag, modelFormat)
        }

    def jobInputs = speciesMetaCh
        .combine( epochInfoCh )
        .map { speciesName, genomeFa, annotGtf,
               trainName, idx, epochDir, evalDir, hmmFlag, modelFormat ->

            tuple(
                trainName, idx,
                epochDir, evalDir, hmmFlag, modelFormat,
                speciesName, genomeFa, annotGtf
            )
        }

    def speciesCsv = speciesList.join(',')
    def mode = splitName  // 'val' or 'test'

    def gffOut = jobInputs | RUN_TIBERIUS | RUN_GFFCOMPARE

    gffOut
        .groupTuple(by: [0, 1])
        .map { trainName, evalDir, statsFiles ->
            tuple(trainName, evalDir, speciesCsv, mode, statsFiles)
        }
        | PLOT_EVAL
}
