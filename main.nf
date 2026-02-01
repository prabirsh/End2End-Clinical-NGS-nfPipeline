#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    Clinical NGS Analysis Pipeline
========================================================================================
    Pipeline Maintainer: Prabir
    
    A production-grade Next Generation Sequencing analysis pipeline for clinical
    applications including germline and somatic variant detection.
    
    Key Features:
    - GATK HaplotypeCaller for germline variant calling
    - GATK Mutect2 for somatic variant calling with Panel of Normals
    - Comprehensive QC and coverage analysis
    - Clinical-grade variant annotation
    - Professional HTML reports
    
    Usage:
        nextflow run main.nf -profile <docker|singularity|conda> \\
            --analysis_type <germline|somatic_paired|somatic_tumor_only> \\
            --sample_sheet samples.csv
========================================================================================
    Maintained by: Prabir
========================================================================================
*/

// Import workflow modules
include { GERMLINE_WORKFLOW } from './workflows/germline'
include { SOMATIC_WORKFLOW } from './workflows/somatic'

// Print pipeline header
def printHeader() {
    log.info """
    ╔═══════════════════════════════════════════════════════════════╗
    ║                                                               ║
    ║           Clinical NGS Analysis Pipeline v2.0                 ║
    ║                                                               ║
    ║              Maintained by: Prabir                            ║
    ║                                                               ║
    ╚═══════════════════════════════════════════════════════════════╝
    
    Pipeline Configuration:
    ─────────────────────────────────────────────────────────────────
    Analysis Type    : ${params.analysis_type}
    Sample Sheet     : ${params.sample_sheet}
    Reference Genome : ${params.reference_genome}
    Target BED       : ${params.target_bed}
    Output Directory : ${params.outdir}
    
    Variant Calling:
    ─────────────────────────────────────────────────────────────────
    Germline Caller  : GATK HaplotypeCaller v4.5+
    Somatic Caller   : GATK Mutect2 v4.5+
    Panel of Normals : ${params.pon ?: 'Not configured (required for tumor-only)'}
    
    ═════════════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow {
    // Print header
    printHeader()
    
    // Validate required parameters
    if (!params.sample_sheet || !params.reference_genome || !params.target_bed) {
        error "ERROR: Required parameters missing. Please provide: --sample_sheet, --reference_genome, --target_bed"
    }
    
    // Parse sample sheet
    // Expected format: sample_id,patient_id,sample_type,fastq_1,fastq_2
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            tuple(
                row.sample_id,
                row.patient_id,
                row.sample_type,
                file(row.fastq_1),
                file(row.fastq_2)
            )
        }
        .set { samples_ch }
    
    // Route to appropriate workflow based on analysis type
    if (params.analysis_type == 'germline') {
        log.info "Starting Germline Workflow with GATK HaplotypeCaller..."
        GERMLINE_WORKFLOW(samples_ch)
        
    } else if (params.analysis_type in ['somatic_paired', 'somatic_tumor_only']) {
        log.info "Starting Somatic Workflow with GATK Mutect2..."
        
        if (params.analysis_type == 'somatic_tumor_only' && !params.pon) {
            error "ERROR: Panel of Normals (--pon) is required for tumor-only analysis"
        }
        
        SOMATIC_WORKFLOW(samples_ch)
        
    } else {
        error "ERROR: Invalid analysis_type. Must be: germline, somatic_paired, or somatic_tumor_only"
    }
}

workflow.onComplete {
    log.info """
    ╔═══════════════════════════════════════════════════════════════╗
    ║                    Pipeline Completed                         ║
    ╚═══════════════════════════════════════════════════════════════╝
    
    Status      : ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
    Duration    : ${workflow.duration}
    Results     : ${params.outdir}
    
    ═════════════════════════════════════════════════════════════════
    Maintained by: Prabir
    """.stripIndent()
}

workflow.onError {
    log.error """
    ╔═══════════════════════════════════════════════════════════════╗
    ║                      Pipeline Error                           ║
    ╚═══════════════════════════════════════════════════════════════╝
    
    Error Message: ${workflow.errorMessage}
    Error Report : ${workflow.errorReport}
    
    ═════════════════════════════════════════════════════════════════
    Maintained by: Prabir
    """.stripIndent()
}
