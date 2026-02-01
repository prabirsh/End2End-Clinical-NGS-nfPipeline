# End to End Clinical NGS Nextflow Pipeline 

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Maintained by: Prabir**
## 🎯 What This Pipeline Does

An automated NGS analysis pipeline that runs end-to-end from FASTQ to results.

- Automated **QC checks** with clear pass/fail customize logic  
- **Standalone HTML report** generation  
- Can be launched via **command line** or **web interface**

### Web Interface (Current Status)
- Web UI is a **template/demo**
- Uses **dummy data** for workflow illustration
- Backend integration and walkthrough video **planned**

---
### ▶️ Pipeline Web Interface Demo

[![Web Interface Demo](https://drive.google.com/thumbnail?id=1ja86ALwu3ySlxKqyYUHFpgUeOXG7BzJO&sz=w1000)](https://drive.google.com/file/d/1ja86ALwu3ySlxKqyYUHFpgUeOXG7BzJO/view)

*Demo video showing the web interface workflow (UI template with dummy data).*



## Featuring:

- ✅ **GATK HaplotypeCaller** for germline variant calling
- ✅ **GATK Mutect2** for somatic variant calling with Panel of Normals support
- ✅ Comprehensive QC and coverage analysis
- ✅ Clinical-grade variant annotation
- ✅ Modern web interface and command-line execution
- ✅ Professional HTML reports

---

## 📊 Pipeline Overview

```
┌────────────────┐
│  FASTQ Files   │  Input: Paired-end reads
└───────┬────────┘
        │
        ▼
┌────────────────┐
│   Raw QC       │  FastQC -> Fail if <10M reads, <80% Q30
└───────┬────────┘  Clinical checkpoint: Reject poor quality
        │
        ▼
┌────────────────┐
│  Alignment     │  BWA-MEM -> SAM -> Sorted BAM
└───────┬────────┘  Mark duplicates, BQSR (if known sites)
        │
        ▼
┌────────────────┐
│  Coverage QC   │  Mean depth, uniformity, dropout detection
└───────┬────────┘  Clinical checkpoint: Ensure ≥100x target coverage
        │
        ├──────────────────┬──────────────────┐
        ▼                  ▼                  ▼
   ┌─────────┐      ┌─────────┐       ┌─────────┐
   │Germline │      │ Somatic │       │ Somatic │
   │         │      │ Paired  │       │Tumor-Only│
   │HaplotypeCaller | Mutect2 │       │Mutect2+PON
   └────┬────┘      └────┬────┘       └────┬────┘
        │                │                  │
        └────────┬───────┴──────────────────┘
                 ▼
        ┌────────────────┐
        │Variant Filtering  VCF -> Filtered VCF
        └───────┬────────┘  Depth, VAF, strand bias, population frequency
                │
                ▼
        ┌────────────────┐
        │  Annotation    │  VEP/ANNOVAR + ClinVar + COSMIC + OncoKB
        └───────┬────────┘
                │
                ▼
        ┌────────────────┐
        │Clinical Report │  HTML/PDF with QC metrics + actionable variants
        └────────────────┘
```

---


## 🚀 Quick Start

## Installation

### A. System Requirements

- **Operating System**: Linux (Ubuntu 20.04+, CentOS 7+) or macOS
- **CPU**: 8+ cores recommended
- **RAM**: 32 GB minimum, 64 GB recommended
- **Storage**: 500 GB+ for analysis workspace

### B. Prerequisites

```bash
# Java (required for Nextflow and GATK)
sudo apt install openjdk-11-jdk

# Python 3.8+
sudo apt install python3 python3-pip

# Docker (recommended)
curl -fsSL https://get.docker.com | sh
```

### C. Automated Installation

```bash
# Download the pipeline
git clone https://github.com/prabirsh/End2End-Clinical-NGS-nfPipeline
cd End2End-Clinical-NGS-nfPipeline

# Run installer
sudo ./install.sh
```

### D. Manual Installation

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Pull Docker image
docker pull broadinstitute/gatk:latest

# Install Python dependencies
pip install flask flask-cors pandas numpy
```

---

### Web Interface
```bash
# Start the web server
cd web
python3 api_server.py

# Open browser to http://localhost:5000
```

### Command Line
```bash
# Germline analysis
nextflow run main.nf \
    --analysis_type germline \
    --sample_sheet samples.csv \
    --reference_genome /path/to/GRCh38.fa \
    --target_bed /path/to/targets.bed

# Somatic paired analysis  
nextflow run main.nf \
    --analysis_type somatic_paired \
    --sample_sheet samples.csv \
    --reference_genome /path/to/GRCh38.fa \
    --target_bed /path/to/targets.bed

# Somatic tumor-only analysis (requires Panel of Normals)
nextflow run main.nf \
    --analysis_type somatic_tumor_only \
    --sample_sheet samples.csv \
    --reference_genome /path/to/GRCh38.fa \
    --target_bed /path/to/targets.bed \
    --pon /path/to/panel_of_normals.vcf.gz
```

---

### 5.2 Navigation

**Dashboard**
- View pipeline statistics
- Monitor active jobs
- See variant calling tool status

**New Analysis**
- Select analysis type
- Upload sample sheet
- Configure parameters
- Start pipeline

**Results**
- Browse completed analyses
- Download VCF files
- View HTML reports

**Logs**
- Real-time progress monitoring
- Error tracking
- Performance metrics

**User Manual**
- In-browser documentation
- Quick reference guide

### 5.3 Running an Analysis

1. Click **"New Analysis"** in sidebar
2. Select analysis type from dropdown
3. Upload or specify sample sheet path
4. Enter reference genome path
5. Enter target BED file path
6. (Tumor-only) Enter Panel of Normals path
7. Click **"Start Pipeline"**
8. Switch to **"Logs"** tab to monitor progress

### 5.4 Viewing Results

1. Navigate to **"Results"** tab
2. Click **"View Report"** for any completed sample
3. Download files using download buttons
4. Review QC metrics and variant calls

---

## 6. Command Line Usage

### 6.1 Basic Syntax

```bash
nextflow run main.nf [OPTIONS]
```

### 6.2 Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--analysis_type` | Type of analysis | `germline`, `somatic_paired`, `somatic_tumor_only` |
| `--sample_sheet` | CSV file with sample info | `samples.csv` |
| `--reference_genome` | Reference FASTA file | `GRCh38.fa` |
| `--target_bed` | Target regions BED file | `targets.bed` |

### 6.3 Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | `results` |
| `--pon` | Panel of Normals (tumor-only) | `null` |
| `--known_sites` | Known variants for BQSR | `null` |
| `--min_coverage` | Minimum coverage threshold | `50` |
| `--min_vaf_germline` | Min VAF for germline | `0.20` |
| `--min_vaf_somatic` | Min VAF for somatic | `0.05` |

---

## Support

For issues, questions, or feature requests:

**Maintained by: Prabir**

---

**Last Updated**: February 2026  
**Version**: 2.0
