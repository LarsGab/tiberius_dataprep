#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Pipeline: GFF3 -> TFRecord split by species
 */

////////////////////////////////////////////////////////////////////////
//                     PARAMETER DEFAULTS                              //
////////////////////////////////////////////////////////////////////////

// Override the whole image with --container docker://my/image:tag if needed.
params.tiberius_version = '2.x'
// Accept the legacy --configYAML flag as an alias for --config_yaml. "First
// assignment wins" in Nextflow 26, so a CLI --config_yaml beats the alias,
// and a CLI --configYAML feeds in via the ?: fallback.
params.config_yaml      = params.containsKey('configYAML') ? params.configYAML : "../config/config.yaml"
// NB: do NOT declare YAML-derived params here. In Nextflow 26 the first
// assignment to a param wins; setting them to null at top level would
// silently shadow the values the workflow loads from the YAML.

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

process GFF3_2_GTF {
    tag { "${species}" }
    container { params.container }
    publishDir { "${params.work_dir}/annot_gtf" }, mode: 'copy', pattern: '*.gtf'

    cpus   1
    memory '4 GB'

    input:
        tuple val(species), val(split), path(annot), path(genome)

    output:
        tuple val(species), val(split),
              path("${species}.gtf"),
              path(genome)

    shell:
    """
    if [[ "${annot}" != *.gtf ]]; then
        # convert non-GTF to GTF
        gffread -T "${annot}" -o "${species}.gtf"
    elif [[ ! -e "${species}.gtf" ]]; then
        # only create the symlink if species.gtf does _not_ yet exist
        ln -s "${annot}" "${species}.gtf"
    fi
    """
}


process REFORMAT_ANNOT {
    tag { "${species}" }
    container { params.container }
    publishDir { "${params.work_dir}/annot_gtf" }, mode: 'copy', pattern: '*.gtf'

    cpus   1
    memory '8 GB'

    input:
        tuple val(species), val(split), path(gtf), path(genome)

    output:
        tuple val(species), val(split),
                path("${species}_reformat.gtf"),
                path(genome)


    script:
    """
    reformat_gtf.py \\
        --input ${gtf} \\
        --out   ${species}_reformat.gtf \\
        --prefix ${species}_
    """
}


process LONGEST_ISOFORM {
    tag { "${species}" }
    container { params.container }
    publishDir { "${params.work_dir}/annot_gtf" }, mode: 'copy', pattern: '*.gtf'

    cpus   1
    memory '4 GB'

    input:
        tuple val(species), val(split), path(gtf), path(genome)

    output:
        tuple val(species), val(split),
              path("${species}_longest.gtf"),
              path(genome)

    script:
    """
    select_single_isoform.py ${gtf} > ${species}_longest.gtf
    """
}


process TFRECORD {
    tag { "${species}" }
    container { params.container }
    storeDir { "cache/TFRECORD/${split}" }
    publishDir { "${params.work_dir}/tfrecords/${split}" }, mode: 'copy', pattern: '*.tfrecords'

    cpus   50
    memory '200 GB'
    time  '24h'

    input:
        tuple val(species), val(split), path(gtf), path(genome)

    output:
        tuple val(species), val(split), path('*.tfrecords')

    script:
    """
    write_tfrecord_species.py \\
        --wsize ${params.chunk_size} \\
        --gtf   ${gtf} \\
        --fasta ${genome} \\
        --out   ${species} \\
        --min_seq_len ${params.min_seq_len}
    """
}


process WRITE_SPECIES_LIST {
    tag { "${split}" }
    container { params.container }
    cpus   1
    memory '1 GB'

    publishDir { "${params.work_dir}/tfrecords/" }, mode: 'copy'

    input:
        tuple val(split), path(cfg_yaml)

    output:
        path "${split}/species.txt"

    script:
    """
    mkdir -p ${split}
    python3 - <<'EOF'
import yaml, pathlib
cfg = yaml.safe_load(open('${cfg_yaml}'))
species = cfg['species_split']['${split}']
species = "" if not species else '\\n'.join(species)
pathlib.Path('${split}/species.txt').write_text(species)
EOF
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

    // Load YAML config. Each assignment is guarded so CLI overrides win.
    def cfg = new groovy.yaml.YamlSlurper().parseText( file(params.config_yaml).text )
    if (!params.work_dir)          params.work_dir    = cfg.work_dir   as String
    if (params.min_seq_len == null) params.min_seq_len = (cfg.min_seq_len ?: 500000) as Integer
    if (params.chunk_size  == null) params.chunk_size  = (cfg.chunk_size  ?: 9999)   as Integer

    // Only needed inside this workflow body, so kept as locals.
    def ANNOT_DIR  = (params.annot_dir  ?: cfg.annot_dir)  as String
    def GENOME_DIR = (params.genome_dir ?: cfg.genome_dir) as String

    // species -> split lookup. Validate that no species appears in more than
    // one split (or more than once within a split).
    def speciesSplits = [:].withDefault { [] }
    ['train','val','test'].each { split ->
        (cfg.species_split[split] ?: []).each { sp -> speciesSplits[sp] << split }
    }
    def duplicates = speciesSplits.findAll { sp, splits -> splits.size() > 1 }
    if (duplicates) {
        def lines = duplicates.collect { sp, splits -> "  ${sp}: ${splits.join(', ')}" }
        error "species_split has species appearing more than once:\n${lines.join('\n')}\nEach species must appear in exactly one split (train, val, or test)."
    }
    def speciesToSplit = speciesSplits.collectEntries { sp, splits -> [(sp): splits[0]] }

    // (split, cfg_yaml) channel for WRITE_SPECIES_LIST
    def cfgCh = Channel.value(file(params.config_yaml))
    def speciesListCh = Channel.of('train','val','test').combine(cfgCh)

    // (species, split, annotFile, genomeFa) channel for the main chain
    def metaCh = Channel.from(speciesToSplit.collect { sp, spl ->
        def annotCandidates = [
            file("${ANNOT_DIR}/${sp}.gtf"),
            file("${ANNOT_DIR}/${sp}.gff3"),
            file("${ANNOT_DIR}/${sp}.gff"),
        ]
        def annotFile = annotCandidates.find { it.exists() }
        if (!annotFile)
            error "No annotation file (.gtf, .gff3 or .gff) found for ${sp} in ${ANNOT_DIR}"

        def genomeCandidates = [
            file("${GENOME_DIR}/${sp}.genome.fa"),
            file("${GENOME_DIR}/${sp}.fa"),
            file("${GENOME_DIR}/${sp}.fasta"),
        ]
        def genomeFile = genomeCandidates.find { it.exists() }
        if (!genomeFile)
            error "No genome FASTA (.genome.fa, .fa or .fasta) found for ${sp} in ${GENOME_DIR}"

        tuple(sp, spl, annotFile, genomeFile)
    })

    speciesListCh | WRITE_SPECIES_LIST

    metaCh                                \
        | GFF3_2_GTF                      \
        | REFORMAT_ANNOT                  \
        | LONGEST_ISOFORM                 \
        | TFRECORD
}
