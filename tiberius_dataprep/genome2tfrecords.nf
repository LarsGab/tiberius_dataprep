#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Pipeline: GFF3 → TFRecord split by species
 * Updated: 2025-05-15
 */

////////////////////////////////////////////////////////////////////////
//                     PARAMETERS & CONFIG                             //
////////////////////////////////////////////////////////////////////////

params.container   = 'docker://larsgabriel23/tiberius:dev'
params.configYAML  = "../config/config_dataprep.yaml"

import groovy.yaml.YamlSlurper

def cfg       = new YamlSlurper().parseText( file(params.configYAML).text )
def ANNOT_DIR = cfg.annot_dir   as String
def GENOME_DIR= cfg.genome_dir  as String
def OUT_DIR   = cfg.work_dir    as String


////////////////////////////////////////////////////////////////////////
//                      BUILD META CHANNEL                            //
////////////////////////////////////////////////////////////////////////

def speciesToSplit = [:]
['train','val','test'].each { split ->
    cfg.species_split[split].each { sp -> speciesToSplit[sp] = split }
}

Channel
    .from( speciesToSplit.collect { sp, spl ->
        tuple(sp, spl,
              file("${ANNOT_DIR}/${sp}.gff3"),
              file("${GENOME_DIR}/${sp}.genome.fa"))
    })
    .set { META_CH }


////////////////////////////////////////////////////////////////////////
//                            PROCESSES                               //
////////////////////////////////////////////////////////////////////////

process GFF3_2_GTF {
    tag "${species}"
    publishDir "${OUT_DIR}/annot_gtf", mode: 'copy', pattern: '*.gtf'  

    cpus   1
    memory '4 GB'

    input:
        tuple val(species), val(split), path(gff3), path(genome)

    output:
        tuple val(species), val(split),
              path("${species}.gtf"),
              path(genome)

    shell:
    """
    gffread -T ${gff3} -o ${species}.gtf
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
    get_longest_isoform.py ${gtf} > ${species}_longest.gtf
    """
}


process TFRECORD {
    tag "${species}"
    container params.container
    publishDir "${OUT_DIR}/tfrecords/${split}", mode: 'move', pattern: '*.tfrecords' 

    cpus   50
    memory '950 GB'

    input:
        tuple val(species), val(split), path(gtf), path(genome)

    output:
        path '*.tfrecords'
        path "${species}_done.txt"

    script:
    """
    write_tfrecord_species.py \
        --wsize 9999 \
        --gtf   ${gtf} \
        --fasta ${genome} \
        --out   ${species} 
    touch ${species}_done.txt
    """
}


process WRITE_SPECIES_LIST {
    tag "${split}"
    cpus   1
    memory '1 GB'

    publishDir "${OUT_DIR}/tfrecords/", mode: 'copy'

    input:
        val split

    output:
        path "${split}/species.txt"

    script:
    """
    mkdir -p ${split}
    python3 - <<'EOF'
import yaml, pathlib
cfg = yaml.safe_load(open('${params.configYAML}'))
species = cfg['species_split']['${split}']
species = "" if not species else '\\n'.join(species)
pathlib.Path('${split}/species.txt').write_text('\\n'.join(species))
EOF
    """
}


////////////////////////////////////////////////////////////////////////
//                            WORKFLOW                                //
////////////////////////////////////////////////////////////////////////

workflow {

    Channel.of('train','val','test') | WRITE_SPECIES_LIST

    META_CH                               \
        | GFF3_2_GTF                      \
        | REFORMAT_ANNOT                  \
        | LONGEST_ISOFORM                 \
        | TFRECORD
}
