#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CNV_COVERAGE {
    tag "$sample_id"
    label 'low_memory'
    publishDir "${params.outdir}/cnv", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path target_bed
    
    output:
    tuple val(sample_id), path("${sample_id}.exon_coverage.txt"), emit: exon_coverage
    
    script:
    """
    # Generate exon-level coverage for CNV calling
    bedtools coverage -a ${target_bed} -b ${bam} -mean > ${sample_id}.exon_coverage.txt
    """
}
