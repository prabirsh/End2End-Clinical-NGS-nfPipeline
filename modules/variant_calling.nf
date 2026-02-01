#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
    Maintained by: Prabir
========================================================================================
    VARIANT CALLING MODULE
========================================================================================
    CLINICAL DECISION POINTS:
    
    GERMLINE vs SOMATIC:
    - Germline: Constitutional variants present in all cells (~50% VAF heterozygous)
    - Somatic: Acquired mutations in tumor cells (VAF depends on purity/clonality)
    
    PAIRED vs TUMOR-ONLY:
    - Paired: Matched normal subtracts germline background (GOLD STANDARD)
    - Tumor-only: Higher false positive rate, requires aggressive filtering
    
    CALLER CHOICE:
    - GATK HaplotypeCaller: Germline variants, local reassembly
    - GATK Mutect2: Somatic variants, tumor-normal or tumor-only
    - FreeBayes: Alternative germline caller, good for indels
========================================================================================
*/

process CALL_VARIANTS_GERMLINE {
    tag "$sample_id"
    label 'variant_calling'
    publishDir "${params.outdir}/vcf/germline", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path target_bed
    
    output:
    tuple val(sample_id), path("${sample_id}.germline.raw.vcf.gz"), path("${sample_id}.germline.raw.vcf.gz.tbi"), emit: vcf
    
    script:
    """
    # CLINICAL CONTEXT: Germline Variant Calling
    # Used for: Hereditary cancer panels, pharmacogenomics, rare disease
    # Expected VAF: ~50% (heterozygous), ~100% (homozygous)
    
    if [ "${params.variant_caller}" = "gatk" ]; then
        # GATK HaplotypeCaller - Industry standard for germline
        gatk HaplotypeCaller \\
            -R ${reference} \\
            -I ${bam} \\
            -L ${target_bed} \\
            -O ${sample_id}.germline.raw.vcf.gz \\
            --min-base-quality-score ${params.min_base_quality} \\
            --standard-min-confidence-threshold-for-calling 10.0 \\
            --dont-use-soft-clipped-bases \\
            --native-pair-hmm-threads ${task.cpus}
        
        echo "✓ Germline variants called with GATK HaplotypeCaller" >&2
        
    elif [ "${params.variant_caller}" = "freebayes" ]; then
        # FreeBayes - Alternative caller, good for complex indels
        freebayes \\
            -f ${reference} \\
            -t ${target_bed} \\
            --min-base-quality ${params.min_base_quality} \\
            --min-mapping-quality ${params.min_mapping_quality} \\
            --min-alternate-fraction ${params.min_vaf_germline} \\
            --min-coverage ${params.min_variant_depth} \\
            ${bam} \\
            | bgzip -c > ${sample_id}.germline.raw.vcf.gz
        
        tabix -p vcf ${sample_id}.germline.raw.vcf.gz
        
        echo "✓ Germline variants called with FreeBayes" >&2
    fi
    
    # Count raw variants
    TOTAL_VARS=\$(zcat ${sample_id}.germline.raw.vcf.gz | grep -v "^#" | wc -l)
    echo "Total raw germline variants: \$TOTAL_VARS" >&2
    
    # CLINICAL EXPECTATION:
    # - WES: 20,000-25,000 SNVs + 2,000-3,000 indels per sample
    # - Targeted panel (50 genes): 50-200 variants
    # - If significantly different, investigate potential issues
    """
}

process CALL_VARIANTS_SOMATIC_PAIRED {
    tag "$patient_id"
    label 'variant_calling'
    publishDir "${params.outdir}/vcf/somatic", mode: 'copy'
    
    input:
    tuple val(patient_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path reference
    path target_bed
    
    output:
    tuple val(patient_id), path("${patient_id}.somatic.raw.vcf.gz"), path("${patient_id}.somatic.raw.vcf.gz.tbi"), emit: vcf
    path "${patient_id}.somatic.stats", emit: stats
    
    script:
    """
    # CLINICAL CONTEXT: Tumor-Normal Paired Somatic Calling
    # GOLD STANDARD for somatic variant detection
    # Normal sample subtracts germline background, reducing false positives
    
    # GATK Mutect2 - Purpose-built for somatic variant calling
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal \$(samtools view -H ${normal_bam} | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | head -1) \\
        -L ${target_bed} \\
        -O ${patient_id}.somatic.raw.vcf.gz \\
        --f1r2-tar-gz ${patient_id}.f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus}
    
    # Learn orientation bias artifacts (FFPE C>T artifacts)
    gatk LearnReadOrientationModel \\
        -I ${patient_id}.f1r2.tar.gz \\
        -O ${patient_id}.read-orientation-model.tar.gz
    
    # Get pileup summaries for contamination estimation
    gatk GetPileupSummaries \\
        -I ${tumor_bam} \\
        -V ${params.gnomad_common ?: params.dbsnp} \\
        -L ${target_bed} \\
        -O ${patient_id}.tumor_pileups.table
    
    gatk GetPileupSummaries \\
        -I ${normal_bam} \\
        -V ${params.gnomad_common ?: params.dbsnp} \\
        -L ${target_bed} \\
        -O ${patient_id}.normal_pileups.table
    
    # Calculate contamination
    gatk CalculateContamination \\
        -I ${patient_id}.tumor_pileups.table \\
        -matched ${patient_id}.normal_pileups.table \\
        -O ${patient_id}.contamination.table \\
        --tumor-segmentation ${patient_id}.segments.table
    
    # Filter Mutect2 calls
    gatk FilterMutectCalls \\
        -R ${reference} \\
        -V ${patient_id}.somatic.raw.vcf.gz \\
        --contamination-table ${patient_id}.contamination.table \\
        --tumor-segmentation ${patient_id}.segments.table \\
        --ob-priors ${patient_id}.read-orientation-model.tar.gz \\
        -O ${patient_id}.somatic.filtered.vcf.gz
    
    # Generate statistics
    bcftools stats ${patient_id}.somatic.filtered.vcf.gz > ${patient_id}.somatic.stats
    
    TOTAL_SOMATIC=\$(zcat ${patient_id}.somatic.filtered.vcf.gz | grep -v "^#" | grep "PASS" | wc -l)
    
    echo "=== Somatic Calling Summary for ${patient_id} ===" >&2
    echo "Total PASS somatic variants: \$TOTAL_SOMATIC" >&2
    echo "" >&2
    echo "CLINICAL EXPECTATION:" >&2
    echo "- Targeted hotspot panel: 1-5 actionable mutations" >&2
    echo "- Comprehensive panel: 5-15 mutations" >&2
    echo "- WES solid tumor: 50-200 mutations (varies by cancer type)" >&2
    echo "- Hypermutated tumors (MSI-H): >1000 mutations" >&2
    """
}

process CALL_VARIANTS_SOMATIC_TUMOR_ONLY {
    tag "$sample_id"
    label 'variant_calling'
    publishDir "${params.outdir}/vcf/somatic", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(sample_type), path(fastq_1), path(fastq_2)
    path reference
    path target_bed
    path pon
    
    output:
    tuple val(sample_id), path("${sample_id}.somatic.raw.vcf.gz"), path("${sample_id}.somatic.raw.vcf.gz.tbi"), emit: vcf
    
    script:
    """
    # CLINICAL CONTEXT: Tumor-Only Somatic Calling
    # Used when matched normal unavailable (FFPE archives, deceased patients)
    # 
    # CRITICAL LIMITATION:
    # Cannot distinguish somatic from germline without matched normal
    # REQUIRES aggressive filtering:
    # - Panel of Normals (PON) from ≥30 normal samples
    # - Population databases (gnomAD, ExAC)
    # - Higher VAF threshold (≥10% vs ≥5% for paired)
    
    echo "⚠️  WARNING: Tumor-only analysis in progress" >&2
    echo "Higher false positive rate expected without matched normal" >&2
    
    # Mutect2 in tumor-only mode
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${bam} \\
        -L ${target_bed} \\
        --panel-of-normals ${pon} \\
        -O ${sample_id}.somatic.raw.vcf.gz \\
        --f1r2-tar-gz ${sample_id}.f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus}
    
    # Learn orientation bias
    gatk LearnReadOrientationModel \\
        -I ${sample_id}.f1r2.tar.gz \\
        -O ${sample_id}.read-orientation-model.tar.gz
    
    # Contamination estimation (tumor-only)
    gatk GetPileupSummaries \\
        -I ${bam} \\
        -V ${params.gnomad_common ?: params.dbsnp} \\
        -L ${target_bed} \\
        -O ${sample_id}.pileups.table
    
    gatk CalculateContamination \\
        -I ${sample_id}.pileups.table \\
        -O ${sample_id}.contamination.table
    
    # Filter with tumor-only specific filters
    gatk FilterMutectCalls \\
        -R ${reference} \\
        -V ${sample_id}.somatic.raw.vcf.gz \\
        --contamination-table ${sample_id}.contamination.table \\
        --ob-priors ${sample_id}.read-orientation-model.tar.gz \\
        -O ${sample_id}.somatic.filtered.vcf.gz
    
    # Additional tumor-only filtering
    bcftools filter \\
        -i 'INFO/AF >= 0.10 && FORMAT/AD[0:1] >= 5' \\
        -O z \\
        -o ${sample_id}.somatic.tumor_only_filtered.vcf.gz \\
        ${sample_id}.somatic.filtered.vcf.gz
    
    tabix -p vcf ${sample_id}.somatic.tumor_only_filtered.vcf.gz
    
    echo "✓ Tumor-only filtering complete" >&2
    echo "Variants with VAF <10% removed (higher threshold for tumor-only)" >&2
    """
}
