#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENERATE_GERMLINE_REPORT {
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(qc_json)
    tuple val(sample_id), path(coverage_summary)
    tuple val(sample_id), path(variants_tsv)
    tuple val(sample_id), path(dup_metrics)
    
    output:
    path "${sample_id}.clinical_report.html", emit: html
    
    script:
    """
    python3 ${projectDir}/scripts/generate_report.py \
        --sample ${sample_id} \
        --qc ${qc_json} \
        --coverage ${coverage_summary} \
        --variants ${variants_tsv} \
        --duplicates ${dup_metrics} \
        --output ${sample_id}.clinical_report.html
    """
}

process GENERATE_SOMATIC_REPORT {
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(qc_json)
    tuple val(sample_id), path(coverage_summary)
    tuple val(sample_id), path(variants_tsv)
    tuple val(sample_id), path(dup_metrics)
    tuple val(sample_id), path(purity)
    
    output:
    path "${sample_id}.clinical_report.html", emit: html
    
    script:
    """
    python3 ${projectDir}/scripts/generate_report.py \
        --sample ${sample_id} \
        --qc ${qc_json} \
        --coverage ${coverage_summary} \
        --variants ${variants_tsv} \
        --duplicates ${dup_metrics} \
        --purity ${purity} \
        --analysis somatic \
        --output ${sample_id}.clinical_report.html
    """
}
