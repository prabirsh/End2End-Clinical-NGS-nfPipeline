#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FILTER_VARIANTS {
    tag "$sample_id"
    label 'low_memory'
    publishDir "${params.outdir}/vcf/filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    val analysis_mode
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), path("${sample_id}.filtered.vcf.gz.tbi"), emit: vcf
    path "${sample_id}.filter_report.txt", emit: report
    
    script:
    if (analysis_mode == 'germline')
        """
        bcftools filter \
            -i 'INFO/DP >= ${params.min_variant_depth} && (GT="0/1" && INFO/AF >= ${params.min_vaf_germline} && INFO/AF <= 0.80) || (GT="1/1" && INFO/AF >= 0.90)' \
            -O z \
            -o ${sample_id}.filtered.vcf.gz \
            ${vcf}
        
        tabix -p vcf ${sample_id}.filtered.vcf.gz
        
        echo "Germline filtering applied: depth >=${params.min_variant_depth}, VAF >=${params.min_vaf_germline}" > ${sample_id}.filter_report.txt
        """
    else
        """
        bcftools filter \
            -i 'INFO/DP >= ${params.min_variant_depth} && INFO/AF >= ${params.min_vaf_somatic} && FORMAT/AD[0:1] >= ${params.min_alt_reads}' \
            -O z \
            -o ${sample_id}.filtered.vcf.gz \
            ${vcf}
        
        tabix -p vcf ${sample_id}.filtered.vcf.gz
        
        echo "Somatic filtering applied: depth >=${params.min_variant_depth}, VAF >=${params.min_vaf_somatic}, alt reads >=${params.min_alt_reads}" > ${sample_id}.filter_report.txt
        """
}
