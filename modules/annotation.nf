#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ANNOTATE_VARIANTS {
    tag "$sample_id"
    label 'annotation'
    publishDir "${params.outdir}/annotated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    val analysis_mode
    
    output:
    tuple val(sample_id), path("${sample_id}.annotated.tsv"), emit: tsv
    path "${sample_id}.annotated.vcf.gz", emit: vcf
    
    script:
    """
    # VEP annotation with clinical databases
    vep \
        --input_file ${vcf} \
        --output_file ${sample_id}.annotated.vcf \
        --format vcf \
        --vcf \
        --everything \
        --assembly ${params.vep_assembly} \
        --species ${params.vep_species} \
        --cache \
        --offline \
        --fork ${task.cpus} \
        --force_overwrite
    
    # Convert to TSV for clinical review
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/DP[\t%GT][\t%AD]\n' \
        ${sample_id}.annotated.vcf > ${sample_id}.annotated.tsv
    
    bgzip ${sample_id}.annotated.vcf
    """
}
