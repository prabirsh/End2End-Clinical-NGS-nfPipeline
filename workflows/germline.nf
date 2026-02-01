#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
    Maintained by: Prabir
========================================================================================
    GERMLINE WORKFLOW
========================================================================================
    Constitutional variant detection workflow for hereditary cancer panels and WES
    
    CLINICAL USE CASES:
    - Lynch syndrome screening
    - BRCA1/2 testing
    - Hereditary cancer panels (50-500 genes)
    - Pharmacogenomics
    
    QC STANDARDS:
    - All target bases must achieve ≥50x coverage
    - Variants called at ≥20% VAF (heterozygous ~50%, homozygous ~100%)
    - Strict filtering against common population variants (gnomAD AF > 1%)
========================================================================================
*/

include { FASTQC } from '../modules/qc'
include { MULTIQC } from '../modules/qc'
include { FASTQC_FILTER } from '../modules/qc'
include { BWA_MEM } from '../modules/alignment'
include { MARK_DUPLICATES } from '../modules/alignment'
include { BQSR } from '../modules/alignment'
include { CALCULATE_COVERAGE } from '../modules/coverage'
include { CALL_VARIANTS_GERMLINE } from '../modules/variant_calling'
include { FILTER_VARIANTS } from '../modules/variant_filtering'
include { ANNOTATE_VARIANTS } from '../modules/annotation'
include { CNV_COVERAGE } from '../modules/cnv'
include { GENERATE_GERMLINE_REPORT } from '../modules/reporting'

workflow GERMLINE_WORKFLOW {
    take:
        samples_ch // tuple(sample_id, patient_id, sample_type, fastq_1, fastq_2)
    
    main:
        // Stage 1: Raw QC
        FASTQC(samples_ch)
        FASTQC_FILTER(FASTQC.out.qc_json)
        
        // Stage 2: Alignment
        BWA_MEM(samples_ch, file(params.reference_genome))
        MARK_DUPLICATES(BWA_MEM.out.bam)
        
        // Stage 3: Base Quality Score Recalibration (BQSR)
        // Clinical rationale: Corrects systematic errors in base quality scores
        // Improves variant calling accuracy, especially for low-frequency variants
        if (params.known_sites) {
            BQSR(
                MARK_DUPLICATES.out.bam,
                file(params.reference_genome),
                file(params.known_sites)
            )
            aligned_bam_ch = BQSR.out.bam
        } else {
            aligned_bam_ch = MARK_DUPLICATES.out.bam
            log.warn "BQSR skipped: no known sites provided. Variant quality may be suboptimal."
        }
        
        // Stage 4: Coverage Analysis
        CALCULATE_COVERAGE(
            aligned_bam_ch,
            file(params.target_bed)
        )
        
        // Stage 5: Germline Variant Calling
        CALL_VARIANTS_GERMLINE(
            aligned_bam_ch,
            file(params.reference_genome),
            file(params.target_bed)
        )
        
        // Stage 6: Variant Filtering
        // Germline-specific filters:
        // - Min depth: 50x (clinical recommendation)
        // - VAF: 0.20-0.80 (heterozygous) or >0.90 (homozygous)
        // - Remove common variants (gnomAD AF > 1%)
        FILTER_VARIANTS(
            CALL_VARIANTS_GERMLINE.out.vcf,
            'germline'
        )
        
        // Stage 7: Annotation
        ANNOTATE_VARIANTS(
            FILTER_VARIANTS.out.vcf,
            'germline'
        )
        
        // Stage 8: CNV-ready coverage matrix
        if (params.enable_cnv_analysis) {
            CNV_COVERAGE(
                aligned_bam_ch,
                file(params.target_bed)
            )
        }
        
        // Stage 9: Generate clinical report
        GENERATE_GERMLINE_REPORT(
            FASTQC.out.qc_json,
            CALCULATE_COVERAGE.out.summary,
            ANNOTATE_VARIANTS.out.tsv,
            MARK_DUPLICATES.out.metrics
        )
    
    emit:
        bam = aligned_bam_ch
        vcf = FILTER_VARIANTS.out.vcf
        annotated = ANNOTATE_VARIANTS.out.tsv
        report = GENERATE_GERMLINE_REPORT.out.html
}
