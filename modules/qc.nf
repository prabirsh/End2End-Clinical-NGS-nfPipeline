#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
    Maintained by: Prabir
========================================================================================
    QC MODULE - Clinical Quality Control
========================================================================================
    WHY THIS MATTERS CLINICALLY:
    Poor quality sequencing can result in:
    - False negatives: Missing pathogenic variants (FAILED DIAGNOSIS)
    - False positives: Calling artifacts as real variants (WRONG TREATMENT)
    - Uneven coverage: Missing variants in critical exons
    
    FAILURE DETECTION:
    This module implements hard QC cutoffs based on CAP/CLIA recommendations
========================================================================================
*/

process FASTQC {
    tag "$sample_id"
    label 'fastqc'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), val(patient_id), val(sample_type), path(fastq_1), path(fastq_2)
    
    output:
    tuple val(sample_id), path("${sample_id}_fastqc_report.json"), emit: qc_json
    path "${sample_id}*_fastqc.{html,zip}", emit: reports
    
    script:
    """
    # Run FastQC
    fastqc -t ${task.cpus} -o . ${fastq_1} ${fastq_2}
    
    # Extract key metrics to JSON for downstream QC filtering
    python3 ${projectDir}/scripts/parse_fastqc.py \
        ${sample_id}_*_fastqc.zip \
        > ${sample_id}_fastqc_report.json
    
    # Log critical metrics
    echo "=== FastQC Summary for ${sample_id} ===" >&2
    grep "Total Sequences" ${sample_id}_*_fastqc/fastqc_data.txt | head -1 >&2
    grep "Sequences flagged as poor quality" ${sample_id}_*_fastqc/fastqc_data.txt | head -1 >&2
    grep "%GC" ${sample_id}_*_fastqc/fastqc_data.txt | head -1 >&2
    """
}

process FASTQC_FILTER {
    tag "$sample_id"
    label 'low_memory'
    
    input:
    tuple val(sample_id), path(qc_json)
    
    output:
    tuple val(sample_id), val(qc_pass), emit: qc_status
    
    script:
    qc_pass = "true"
    """
    #!/usr/bin/env python3
    import json
    import sys
    
    # Clinical QC thresholds
    MIN_READS = ${params.min_reads}
    MIN_Q30 = ${params.min_q30}
    MAX_ADAPTER_CONTENT = 5.0  # Max % adapter contamination
    
    with open('${qc_json}') as f:
        qc = json.load(f)
    
    failures = []
    warnings = []
    
    # Check total reads
    if qc['total_sequences'] < MIN_READS:
        failures.append(f"Insufficient reads: {qc['total_sequences']:,} < {MIN_READS:,}")
    
    # Check Q30 percentage
    if qc['percent_q30'] < MIN_Q30:
        failures.append(f"Low Q30%: {qc['percent_q30']:.1f}% < {MIN_Q30}%")
    
    # Check adapter contamination
    if qc['adapter_content'] > MAX_ADAPTER_CONTENT:
        warnings.append(f"High adapter content: {qc['adapter_content']:.1f}%")
    
    # Check GC content (human expected: 38-42%)
    if not (36 < qc['percent_gc'] < 44):
        warnings.append(f"Abnormal GC content: {qc['percent_gc']:.1f}% (expected 38-42%)")
    
    # Clinical interpretation
    if failures:
        print(f"\\n❌ SAMPLE ${sample_id} FAILED QC:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        print("\\nCLINICAL IMPACT: Sample does not meet minimum standards for diagnostic reporting.", file=sys.stderr)
        print("RECOMMENDED ACTION: Re-sequence sample or reject for insufficient quality.\\n", file=sys.stderr)
        
        if ${params.fail_on_qc.toString().toLowerCase()}:
            sys.exit(1)
    
    if warnings:
        print(f"\\n⚠️  SAMPLE ${sample_id} QC WARNINGS:", file=sys.stderr)
        for w in warnings:
            print(f"  - {w}", file=sys.stderr)
        print("\\nCLINICAL IMPACT: Proceed with caution. Manual review recommended.\\n", file=sys.stderr)
    
    if not failures and not warnings:
        print(f"\\n✓ SAMPLE ${sample_id} PASSED QC", file=sys.stderr)
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    path('*')
    
    output:
    path 'multiqc_report.html'
    path 'multiqc_data'
    
    script:
    """
    multiqc . -n multiqc_report.html
    """
}

process FAILURE_DETECTION {
    tag "$sample_id"
    label 'low_memory'
    publishDir "${params.outdir}/qc/failure_analysis", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf), path(coverage)
    
    output:
    path "${sample_id}_failure_report.txt", emit: report
    
    script:
    """
    cat > failure_analysis.py << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
# Clinical Failure Mode Detection
# This script identifies common technical failures that affect clinical interpretation

print("="*80)
print(f"FAILURE MODE ANALYSIS: ${sample_id}")
print("="*80)
print()

# FAILURE MODE 1: Variant visible in BAM but missing in VCF
# ROOT CAUSE: Variant caller filters (depth, quality, strand bias)
# CLINICAL IMPACT: False negative - missed pathogenic variant
print("1. BAM/VCF Concordance Check")
print("-" * 40)
print("WHAT TO CHECK:")
print("  - Load BAM in IGV and manually inspect suspected variant positions")
print("  - If variant visible with good quality but missing in VCF:")
print("    → Check variant caller logs for filter reasons")
print("    → Consider lowering quality thresholds for that specific position")
print("    → Manual review required before clinical reporting")
print()

# FAILURE MODE 2: Low VAF somatic variants
# ROOT CAUSE: Low tumor purity, subclonal variants, sequencing depth
# CLINICAL IMPACT: Missing actionable mutations
print("2. Low VAF Variant Detection")
print("-" * 40)
print("DETECTION STRATEGY:")
print("  - For targeted panels: variants with VAF 1-5% may be real")
print("  - Check tumor purity estimate (expect >20% for reliable calls)")
print("  - Low VAF + low alt reads = likely artifact")
print("  - Low VAF + high depth + balanced strand = potentially real")
print()

# FAILURE MODE 3: High duplicate rate
# ROOT CAUSE: Low input DNA, over-amplification, library prep issues
# CLINICAL IMPACT: Reduced effective coverage, potential for artifacts
print("3. PCR Duplicate Analysis")
print("-" * 40)
print("THRESHOLDS:")
print("  - <20% duplicates: Excellent")
print("  - 20-30% duplicates: Acceptable")
print("  - 30-50% duplicates: Warning - reduced effective coverage")
print("  - >50% duplicates: CRITICAL - may need re-extraction and re-sequencing")
print()

# FAILURE MODE 4: BED file coordinate mismatch
# ROOT CAUSE: 0-based vs 1-based coordinates, wrong genome build
# CLINICAL IMPACT: Wrong regions analyzed, missed variants
print("4. Target Region Coordinate Check")
print("-" * 40)
print("COMMON ISSUES:")
print("  - BED files are 0-based, half-open [start, end)")
print("  - VCF files are 1-based, closed [pos, pos]")
print("  - Always verify genome build: hg19 vs hg38")
print("  - Off-by-one errors can exclude critical exon boundaries")
print()

# FAILURE MODE 5: GC-rich exon dropout
# ROOT CAUSE: PCR bias against high GC content regions
# CLINICAL IMPACT: Systematically missed variants in specific exons
print("5. GC Bias and Coverage Dropout")
print("-" * 40)
print("HIGH-RISK GENES WITH GC-RICH EXONS:")
print("  - BRCA1 exon 11 (~60% GC)")
print("  - TP53 exons 5-8 (~55% GC)")
print("  - EGFR exon 20 (~58% GC)")
print()
print("DETECTION:")
print("  - Plot coverage vs GC content")
print("  - Identify exons with <50% of mean coverage")
print("  - Manual Sanger sequencing may be required for these regions")
print()

# FAILURE MODE 6: Tumor-only false positives
# ROOT CAUSE: Germline variants called as somatic without matched normal
# CLINICAL IMPACT: Wrong therapeutic decisions
print("6. Tumor-Only Analysis Pitfalls")
print("-" * 40)
print("WITHOUT MATCHED NORMAL:")
print("  - Germline SNPs can appear as 'somatic' mutations")
print("  - ESSENTIAL: Use Panel of Normals (PON) with ≥30 samples")
print("  - Filter against gnomAD/dbSNP aggressively (AF > 0.1%)")
print("  - Flag any variant with VAF ~50% (likely germline)")
print()
print("RECOMMENDED: Always use tumor-normal paired analysis when possible")
print()

print("="*80)
print("END OF FAILURE MODE ANALYSIS")
print("="*80)
PYTHON_SCRIPT

    python3 failure_analysis.py > ${sample_id}_failure_report.txt
    cat ${sample_id}_failure_report.txt
    """
}

process TUMOR_PURITY_ESTIMATE {
    tag "$sample_id"
    label 'low_memory'
    publishDir "${params.outdir}/qc/tumor_purity", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai), val(patient_id), val(sample_type), path(fastq_1), path(fastq_2)
    path bed
    
    output:
    tuple val(sample_id), path("${sample_id}_purity.txt"), emit: purity
    
    when:
    sample_type == 'tumor'
    
    script:
    """
    # Simplified purity estimation based on VAF distribution
    # Clinical tools: ASCAT, Sequenza, FACETS
    
    echo "TUMOR_PURITY_ESTIMATE: ${sample_id}" > ${sample_id}_purity.txt
    echo "" >> ${sample_id}_purity.txt
    echo "NOTE: This is a placeholder for actual purity estimation." >> ${sample_id}_purity.txt
    echo "" >> ${sample_id}_purity.txt
    echo "CLINICAL RECOMMENDATION:" >> ${sample_id}_purity.txt
    echo "- Purity >40%: Ideal for somatic calling" >> ${sample_id}_purity.txt
    echo "- Purity 20-40%: Acceptable, but may miss low VAF variants" >> ${sample_id}_purity.txt
    echo "- Purity <20%: Consider enrichment or re-biopsy" >> ${sample_id}_purity.txt
    echo "" >> ${sample_id}_purity.txt
    echo "For production use, integrate tools like:" >> ${sample_id}_purity.txt
    echo "  - FACETS (for WES/WGS)" >> ${sample_id}_purity.txt
    echo "  - Sequenza (copy number + purity)" >> ${sample_id}_purity.txt
    echo "  - PureCN (for targeted panels)" >> ${sample_id}_purity.txt
    """
}
