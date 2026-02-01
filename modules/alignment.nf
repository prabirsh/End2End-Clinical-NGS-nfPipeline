#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
    Maintained by: Prabir
========================================================================================
    ALIGNMENT MODULE
========================================================================================
    CLINICAL RATIONALE:
    Accurate read alignment is foundational to variant detection. Misalignment can:
    - Create false positive variants (reads mapped to wrong location)
    - Miss real variants (reads unmapped or misaligned)
    - Affect CNV analysis (incorrect depth calculations)
    
    TOOL CHOICE: BWA-MEM
    - Industry standard for clinical NGS
    - DRAGEN-compatible (can swap to DRAGEN in production)
    - Handles both short reads (WES/panels) and longer reads (NovaSeq)
========================================================================================
*/

process BWA_MEM {
    tag "$sample_id"
    label 'alignment'
    publishDir "${params.outdir}/bam", mode: 'copy', enabled: params.keep_intermediate
    
    input:
    tuple val(sample_id), val(patient_id), val(sample_type), path(fastq_1), path(fastq_2)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.alignment_metrics.txt"), emit: metrics
    
    script:
    // Read group information - CRITICAL for clinical pipelines
    // Required by GATK and for merging multiple lanes/flowcells
    def rg_id = sample_id
    def rg_sm = sample_id
    def rg_lb = "lib_${sample_id}"
    def rg_pl = "ILLUMINA"
    def rg_pu = "unit1"
    
    """
    # Clinical alignment pipeline with proper read group annotation
    
    # WHY READ GROUPS MATTER:
    # - Sample tracking across merge operations
    # - Quality score recalibration (BQSR) requires RG info
    # - Duplicate marking across lanes
    # - Clinical QC and chain of custody
    
    bwa mem \\
        -t ${task.cpus} \\
        -M \\
        -R '@RG\\tID:${rg_id}\\tSM:${rg_sm}\\tLB:${rg_lb}\\tPL:${rg_pl}\\tPU:${rg_pu}' \\
        ${reference} \\
        ${fastq_1} \\
        ${fastq_2} \\
        | samtools sort -@ ${task.cpus} -m 2G -o ${sample_id}.sorted.bam -
    
    # Index immediately
    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    
    # Generate alignment metrics
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.alignment_metrics.txt
    samtools stats ${sample_id}.sorted.bam >> ${sample_id}.alignment_metrics.txt
    
    # Log critical metrics
    echo "=== Alignment Summary for ${sample_id} ===" >&2
    grep "mapped (" ${sample_id}.alignment_metrics.txt | head -1 >&2
    grep "properly paired" ${sample_id}.alignment_metrics.txt | head -1 >&2
    
    # CLINICAL QC CHECKPOINT
    MAPPED_READS=\$(grep "mapped (" ${sample_id}.alignment_metrics.txt | head -1 | cut -f1 -d' ')
    TOTAL_READS=\$(grep "in total" ${sample_id}.alignment_metrics.txt | head -1 | cut -f1 -d' ')
    MAPPING_RATE=\$(echo "scale=2; \$MAPPED_READS * 100 / \$TOTAL_READS" | bc)
    
    echo "Mapping rate: \$MAPPING_RATE%" >&2
    
    # Fail if mapping rate too low
    if (( \$(echo "\$MAPPING_RATE < 80" | bc -l) )); then
        echo "ERROR: Mapping rate \$MAPPING_RATE% < 80%" >&2
        echo "CLINICAL IMPACT: Low mapping rate indicates:" >&2
        echo "  - Wrong reference genome (hg19 vs hg38)" >&2
        echo "  - Poor sample quality" >&2
        echo "  - Contamination with non-human DNA" >&2
        echo "  - Adapter sequences not trimmed" >&2
        exit 1
    fi
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    label 'alignment'
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.dedup_metrics.txt"), emit: metrics
    
    script:
    """
    # CLINICAL RATIONALE FOR DUPLICATE MARKING:
    # Optical/PCR duplicates artificially inflate coverage and can:
    # - Create false positive variants (same error amplified)
    # - Skew VAF calculations (artifacts appear at higher frequency)
    # - Violate statistical assumptions in variant callers
    
    # Mark duplicates (Picard-compatible)
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.dedup.bam \\
        -M ${sample_id}.dedup_metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 \\
        --ASSUME_SORT_ORDER coordinate
    
    # Parse duplicate metrics
    DUP_RATE=\$(grep -A 1 "^LIBRARY" ${sample_id}.dedup_metrics.txt | tail -1 | cut -f9)
    
    echo "=== Duplicate Metrics for ${sample_id} ===" >&2
    echo "Duplicate rate: \$DUP_RATE" >&2
    
    # CLINICAL INTERPRETATION
    if (( \$(echo "\$DUP_RATE > 0.50" | bc -l) )); then
        echo "⚠️  CRITICAL: Duplicate rate > 50%" >&2
        echo "CAUSE: Likely low input DNA or excessive PCR cycles" >&2
        echo "IMPACT: Effective coverage significantly reduced" >&2
        echo "ACTION: Consider re-extraction and library prep" >&2
    elif (( \$(echo "\$DUP_RATE > 0.30" | bc -l) )); then
        echo "⚠️  WARNING: Duplicate rate > 30%" >&2
        echo "IMPACT: Moderate reduction in effective coverage" >&2
        echo "ACTION: Review library prep protocol" >&2
    else
        echo "✓ Acceptable duplicate rate" >&2
    fi
    """
}

process BQSR {
    tag "$sample_id"
    label 'alignment'
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.recal.bam*"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path known_sites
    
    output:
    tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bam.bai"), emit: bam
    path "${sample_id}.recal_data.table", emit: recal_table
    
    script:
    """
    # CLINICAL RATIONALE FOR BASE QUALITY SCORE RECALIBRATION (BQSR):
    # 
    # Sequencers systematically over/under-estimate base quality scores based on:
    # - Machine cycle (quality degrades toward end of read)
    # - Sequence context (homopolymers have higher error)
    # - Dinucleotide context (GG has different error than AT)
    # 
    # WHY THIS MATTERS:
    # - Variant callers rely on quality scores for confidence
    # - Uncalibrated scores lead to false positives AND false negatives
    # - Critical for low VAF somatic variant detection
    # 
    # LIMITATION:
    # - Requires known variant sites (dbSNP, 1000G, Mills indels)
    # - Cannot correct for novel systematic errors
    
    # Step 1: Build recalibration model
    gatk BaseRecalibrator \\
        -I ${bam} \\
        -R ${reference} \\
        --known-sites ${known_sites} \\
        -O ${sample_id}.recal_data.table
    
    # Step 2: Apply recalibration
    gatk ApplyBQSR \\
        -I ${bam} \\
        -R ${reference} \\
        --bqsr-recal-file ${sample_id}.recal_data.table \\
        -O ${sample_id}.recal.bam \\
        --create-output-bam-index true
    
    echo "✓ BQSR completed for ${sample_id}" >&2
    echo "Quality scores recalibrated based on known variant sites" >&2
    """
}
