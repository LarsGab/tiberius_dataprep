#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Pipeline: GFF3 → TFRecord split by species
 * Updated: 2025-05-15
 */

////////////////////////////////////////////////////////////////////////
//                     PARAMETERS & CONFIG                             //
////////////////////////////////////////////////////////////////////////

// Supported Tiberius versions -> pinned container tag.
// Override the whole image with --container docker://my/image:tag if needed.
def TIBERIUS_TAG_FOR = [ '1.x': '1.1.8', '2.x': '2.0.6' ]
params.tiberius_version = '2.x'
if (!TIBERIUS_TAG_FOR.containsKey(params.tiberius_version)) {
    error "params.tiberius_version must be one of ${TIBERIUS_TAG_FOR.keySet()} (got '${params.tiberius_version}')"
}
params.container   = "docker://larsgabriel23/tiberius:${TIBERIUS_TAG_FOR[params.tiberius_version]}"
params.config_yaml  = "../config/config.yaml"

import groovy.yaml.YamlSlurper

def cfg       = new YamlSlurper().parseText( file(params.config_yaml).text )
def ANNOT_DIR = cfg.annot_dir   as String
def GENOME_DIR= cfg.genome_dir  as String
def OUT_DIR   = cfg.work_dir    as String
// Minimum genome sequence length (in bp) eligible for training.
// Sequences shorter than this are filtered out before being chunked.
def MIN_SEQ_LEN = cfg.min_seq_len ?: 500000
// Size (in bp) of each TFRecord chunk written from a genome sequence.
// Forwarded to write_tfrecord_species.py as --wsize.
def CHUNK_SIZE  = cfg.chunk_size  ?: 9999
def CFG_CH = Channel.value( file(params.config_yaml) )

////////////////////////////////////////////////////////////////////////
//                      BUILD META CHANNEL                            //
////////////////////////////////////////////////////////////////////////

def speciesToSplit = [:]
['train','val','test'].each { split ->
    cfg.species_split[split].each { sp -> speciesToSplit[sp] = split }
}

Channel
    .of('train','val','test')
    .combine(CFG_CH)      
    .set   { SPECIES_LIST_CH }         // (split , cfg_yaml)

Channel
    .from( speciesToSplit.collect { sp, spl ->
        def gtfFile  = file("${ANNOT_DIR}/${sp}.gtf")
        def gff3File = file("${ANNOT_DIR}/${sp}.gff3")
        def gffFile = file("${ANNOT_DIR}/${sp}.gff")

        def annotFile
        if (gtfFile.exists()) {
            annotFile = gtfFile
        } else if (gff3File.exists()) {
            annotFile = gff3File
        } else if (gffFile.exists()) {
            annotFile = gffFile
        } else {
            error "No annotation file (.gtf, .gff or .gff3) found for ${sp} in ${ANNOT_DIR}"
        }

        tuple( sp, spl, annotFile,
            file("${GENOME_DIR}/${sp}.genome.fa") )
    })
    .set { META_CH }
////////////////////////////////////////////////////////////////////////
//                            PROCESSES                               //
////////////////////////////////////////////////////////////////////////

process GFF3_2_GTF {
    tag "${species}"
    container params.container
    publishDir "${OUT_DIR}/annot_gtf", mode: 'copy', pattern: '*.gtf'

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
    tag "${species}"
    container params.container
    publishDir "${OUT_DIR}/annot_gtf", mode: 'copy', pattern: '*.gtf'  

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
    reformat_gtf.py \
        --input ${gtf} \
        --out   ${species}_reformat.gtf \
        --prefix ${species}_
    """
}


process LONGEST_ISOFORM {
    tag "${species}"
    container params.container
    publishDir "${OUT_DIR}/annot_gtf", mode: 'copy', pattern: '*.gtf'  

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
    tag "${species}"
    container params.container
    storeDir "cache/TFRECORD/${split}"
    publishDir "${OUT_DIR}/tfrecords/${split}", mode: 'copy', pattern: '*.tfrecords'

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
        --wsize ${CHUNK_SIZE} \\
        --gtf   ${gtf} \\
        --fasta ${genome} \\
        --out   ${species} \\
        --min_seq_len ${MIN_SEQ_LEN}
    """
}


process WRITE_SPECIES_LIST {
    tag "${split}"
    container params.container
    cpus   1
    memory '1 GB'

    publishDir "${OUT_DIR}/tfrecords/", mode: 'copy'

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

    SPECIES_LIST_CH  | WRITE_SPECIES_LIST

    META_CH                               \
        | GFF3_2_GTF                      \
        | REFORMAT_ANNOT                  \
        | LONGEST_ISOFORM                 \
        | TFRECORD
}
