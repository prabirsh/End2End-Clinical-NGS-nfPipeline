#!/bin/bash
#
# Run Pipeline from Terminal
# Maintained by: Prabir
#
# Example script for running the pipeline from command line

set -e

echo "╔════════════════════════════════════════════════════════╗"
echo "║                                                        ║"
echo "║      Clinical NGS Pipeline - Terminal Mode            ║"
echo "║              Maintained by: Prabir                     ║"
echo "║                                                        ║"
echo "╚════════════════════════════════════════════════════════╝"
echo ""

# Default values
ANALYSIS_TYPE="germline"
SAMPLE_SHEET=""
REFERENCE=""
TARGET_BED=""
PON=""
PROFILE="docker"
OUTDIR="results"

# Print usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Required Options:"
    echo "  -t, --type TYPE          Analysis type: germline, somatic_paired, somatic_tumor_only"
    echo "  -s, --samples FILE       Sample sheet CSV file"
    echo "  -r, --reference FILE     Reference genome FASTA file"
    echo "  -b, --bed FILE           Target regions BED file"
    echo ""
    echo "Optional:"
    echo "  -p, --pon FILE           Panel of Normals (required for tumor-only)"
    echo "  -o, --outdir DIR         Output directory (default: results)"
    echo "  --profile PROFILE        Execution profile: docker, singularity, conda (default: docker)"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Germline analysis"
    echo "  $0 -t germline -s samples.csv -r ref.fa -b targets.bed"
    echo ""
    echo "  # Somatic paired"
    echo "  $0 -t somatic_paired -s tumor_normal.csv -r ref.fa -b targets.bed"
    echo ""
    echo "  # Somatic tumor-only"
    echo "  $0 -t somatic_tumor_only -s tumor.csv -r ref.fa -b targets.bed -p pon.vcf.gz"
    echo ""
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--type)
            ANALYSIS_TYPE="$2"
            shift 2
            ;;
        -s|--samples)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -b|--bed)
            TARGET_BED="$2"
            shift 2
            ;;
        -p|--pon)
            PON="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$SAMPLE_SHEET" ] || [ -z "$REFERENCE" ] || [ -z "$TARGET_BED" ]; then
    echo "❌ Error: Missing required arguments"
    echo ""
    usage
fi

# Validate analysis type
if [[ ! "$ANALYSIS_TYPE" =~ ^(germline|somatic_paired|somatic_tumor_only)$ ]]; then
    echo "❌ Error: Invalid analysis type: $ANALYSIS_TYPE"
    echo "   Must be: germline, somatic_paired, or somatic_tumor_only"
    exit 1
fi

# Check for PoN if tumor-only
if [ "$ANALYSIS_TYPE" == "somatic_tumor_only" ] && [ -z "$PON" ]; then
    echo "❌ Error: Panel of Normals (--pon) is required for tumor-only analysis"
    exit 1
fi

# Print configuration
echo "Configuration:"
echo "─────────────────────────────────────────────────────────"
echo "Analysis Type   : $ANALYSIS_TYPE"
echo "Sample Sheet    : $SAMPLE_SHEET"
echo "Reference       : $REFERENCE"
echo "Target BED      : $TARGET_BED"
if [ -n "$PON" ]; then
    echo "Panel of Normals: $PON"
fi
echo "Output Directory: $OUTDIR"
echo "Profile         : $PROFILE"
echo "─────────────────────────────────────────────────────────"
echo ""

# Build command
CMD="nextflow run main.nf \
    -profile $PROFILE \
    --analysis_type $ANALYSIS_TYPE \
    --sample_sheet $SAMPLE_SHEET \
    --reference_genome $REFERENCE \
    --target_bed $TARGET_BED \
    --outdir $OUTDIR"

if [ -n "$PON" ]; then
    CMD="$CMD --pon $PON"
fi

# Confirm before running
read -p "🚀 Start pipeline? [Y/n] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]] && [[ -n $REPLY ]]; then
    echo "Cancelled."
    exit 0
fi

echo ""
echo "Starting pipeline..."
echo "═══════════════════════════════════════════════════════════"
echo ""

# Run pipeline
eval $CMD

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "✅ Pipeline completed!"
echo "Results available in: $OUTDIR"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Maintained by: Prabir"
