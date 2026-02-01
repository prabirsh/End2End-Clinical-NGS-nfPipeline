#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
    Maintained by: Prabir
========================================================================================
    SOMATIC WORKFLOW
========================================================================================
    Tumor-specific variant detection for actionable mutations
    
    CLINICAL USE CASES:
    - Lung/colorectal/breast cancer hotspot panels
    - Comprehensive solid tumor panels (300-500 genes)
    - Liquid biopsy (ctDNA) with ultra-low VAF detection
    
    CRITICAL DIFFERENCES FROM GERMLINE:
    - Lower VAF threshold (1-5% for targeted panels, 10% for tumor-only)
    - Tumor-normal paired analysis preferred
    - Panel of Normals (PON) filtering for tumor-only
    - Strand bias and context-specific error filtering
    
    FAILURE MODES TO WATCH:
    - Tumor purity <20%: Low VAF variants may be missed
    - High tumor necrosis: Coverage dropout in specific regions
    - FFPE artifacts: C>T transitions require aggressive filtering
========================================================================================
*/

include { FASTQC } from '../modules/qc'
include { MULTIQC } from '../modules/qc'
include { FASTQC_FILTER } from '../modules/qc'
include { BWA_MEM } from '../modules/alignment'
include { MARK_DUPLICATES } from '../modules/alignment'
include { CALCULATE_COVERAGE } from '../modules/coverage'
include { CALL_VARIANTS_SOMATIC_PAIRED } from '../modules/variant_calling'
include { CALL_VARIANTS_SOMATIC_TUMOR_ONLY } from '../modules/variant_calling'
include { FILTER_VARIANTS } from '../modules/variant_filtering'
include { ANNOTATE_VARIANTS } from '../modules/annotation'
include { CNV_COVERAGE } from '../modules/cnv'
include { GENERATE_SOMATIC_REPORT } from '../modules/reporting'
include { TUMOR_PURITY_ESTIMATE } from '../modules/qc'

workflow SOMATIC_WORKFLOW {
    take:
        samples_ch // tuple(sample_id, patient_id, sample_type, fastq_1, fastq_2)
    
    main:
        // Stage 1: Raw QC
        FASTQC(samples_ch)
        FASTQC_FILTER(FASTQC.out.qc_json)
        
        // Stage 2: Alignment
        BWA_MEM(samples_ch, file(params.reference_genome))
        MARK_DUPLICATES(BWA_MEM.out.bam)
        
        // Stage 3: Coverage Analysis
        CALCULATE_COVERAGE(
            MARK_DUPLICATES.out.bam,
            file(params.target_bed)
        )
        
        // Stage 4: Tumor Purity Estimation
        // Clinical rationale: Low purity samples may need adjusted VAF thresholds
        tumor_samples = samples_ch.filter { it[2] == 'tumor' }
        TUMOR_PURITY_ESTIMATE(
            MARK_DUPLICATES.out.bam.join(tumor_samples),
            file(params.target_bed)
        )
        
        // Stage 5: Somatic Variant Calling
        if (params.analysis_type == 'somatic_paired') {
            // Tumor-Normal Paired Analysis (GOLD STANDARD)
            // Group by patient_id to pair tumor and normal
            paired_ch = MARK_DUPLICATES.out.bam
                .map { sample_id, bam, bai ->
                    def patient_id = sample_id.split('_')[0] // assumes PATIENT_tumor/PATIENT_normal naming
                    tuple(patient_id, sample_id, bam, bai)
                }
                .groupTuple(by: 0)
                .map { patient_id, sample_ids, bams, bais ->
                    def tumor_idx = sample_ids.findIndexOf { it.contains('tumor') }
                    def normal_idx = sample_ids.findIndexOf { it.contains('normal') }
                    
                    if (tumor_idx == -1 || normal_idx == -1) {
                        log.error "Missing tumor or normal for patient: ${patient_id}"
                        return null
                    }
                    
                    tuple(
                        patient_id,
                        bams[tumor_idx], bais[tumor_idx],
                        bams[normal_idx], bais[normal_idx]
                    )
                }
                .filter { it != null }
            
            CALL_VARIANTS_SOMATIC_PAIRED(
                paired_ch,
                file(params.reference_genome),
                file(params.target_bed)
            )
            
            variants_ch = CALL_VARIANTS_SOMATIC_PAIRED.out.vcf
            
        } else {
            // Tumor-Only Analysis (requires Panel of Normals)
            // Clinical context: Used when matched normal unavailable
            // CRITICAL: Higher false positive rate without matched normal
            
            if (!params.pon) {
                log.error "ERROR: Tumor-only analysis requires Panel of Normals (--pon)"
                System.exit(1)
            }
            
            tumor_bams = MARK_DUPLICATES.out.bam.join(tumor_samples)
            
            CALL_VARIANTS_SOMATIC_TUMOR_ONLY(
                tumor_bams,
                file(params.reference_genome),
                file(params.target_bed),
                file(params.pon)
            )
            
            variants_ch = CALL_VARIANTS_SOMATIC_TUMOR_ONLY.out.vcf
        }
        
        // Stage 6: Variant Filtering
        // Somatic-specific filters:
        // - Min VAF: 5% (paired) or 10% (tumor-only)
        // - Min alt reads: 5
        // - Strand bias filter
        // - Clustered variant filter (FFPE artifact detection)
        FILTER_VARIANTS(
            variants_ch,
            'somatic'
        )
        
        // Stage 7: Annotation
        // Prioritize:
        // - COSMIC mutations (known oncogenic)
        // - ClinVar pathogenic variants
        // - OncoKB actionable mutations
        ANNOTATE_VARIANTS(
            FILTER_VARIANTS.out.vcf,
            'somatic'
        )
        
        // Stage 8: CNV Analysis
        if (params.enable_cnv_analysis) {
            CNV_COVERAGE(
                MARK_DUPLICATES.out.bam,
                file(params.target_bed)
            )
        }
        
        // Stage 9: Generate Clinical Report
        GENERATE_SOMATIC_REPORT(
            FASTQC.out.qc_json,
            CALCULATE_COVERAGE.out.summary,
            ANNOTATE_VARIANTS.out.tsv,
            MARK_DUPLICATES.out.metrics,
            TUMOR_PURITY_ESTIMATE.out.purity
        )
    
    emit:
        bam = MARK_DUPLICATES.out.bam
        vcf = FILTER_VARIANTS.out.vcf
        annotated = ANNOTATE_VARIANTS.out.tsv
        report = GENERATE_SOMATIC_REPORT.out.html
}
