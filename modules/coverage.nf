#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CALCULATE_COVERAGE {
    tag "$sample_id"
    label 'low_memory'
    publishDir "${params.outdir}/coverage", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path target_bed
    
    output:
    tuple val(sample_id), path("${sample_id}.coverage_summary.txt"), emit: summary
    path "${sample_id}.coverage_per_base.txt", emit: per_base
    
    script:
    """
    # Calculate coverage metrics
    bedtools coverage -a ${target_bed} -b ${bam} -d > ${sample_id}.coverage_per_base.txt
    
    # Summary statistics
    echo "=== Coverage Summary for ${sample_id} ===" > ${sample_id}.coverage_summary.txt
    samtools depth -b ${target_bed} ${bam} | \
        awk '{sum+=\$3; if(\$3>=100) count++} END {print "Mean depth: "sum/NR; print "Bases >=100x: "count}' \
        >> ${sample_id}.coverage_summary.txt
    """
}
