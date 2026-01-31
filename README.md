# Clinical NGS Pipeline - Production-Grade Oncology Diagnostics

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🎯 What Problem This Solves

**Clinical Reality**: Oncology diagnostic labs process hundreds of NGS samples monthly for actionable mutation detection. Each sample must meet strict CAP/CLIA quality standards, and every variant call directly impacts treatment decisions (targeted therapies, immunotherapy eligibility, clinical trial enrollment).

**The Challenge**: 
- Manual pipelines are error-prone and non-reproducible
- Commercial platforms (Foundation Medicine, Guardant) are expensive (~$3000/test)
- Academic pipelines lack clinical rigor (no QC enforcement, no failure detection)

**This Pipeline Solves**:
✅ End-to-end automation from FASTQ → clinically-actionable variant report  
✅ Enforced QC thresholds aligned with CAP/CLIA guidelines  
✅ Intelligent failure detection (not just error codes, but *why* it failed)  
✅ Somatic + germline workflows with proper filtering strategies  
✅ Production-ready: scalable to 100s of samples, cloud-compatible  

---

## 🏥 Clinical Use Cases

### Primary Applications
1. **Hereditary Cancer Testing**
   - BRCA1/2, Lynch syndrome, Li-Fraumeni
   - Requires: 50x coverage, strict germline filtering
   
2. **Solid Tumor Profiling**
   - Lung (EGFR, ALK, ROS1), colorectal (KRAS, NRAS, BRAF)
   - Requires: Tumor-normal pairing or robust PON

3. **Liquid Biopsy (ctDNA)**
   - Ultra-low VAF detection (1-5%)
   - Requires: High depth (500-1000x), error suppression

4. **Pharmacogenomics**
   - DPYD, TPMT, UGT1A1 variants for chemotherapy dosing
   - Requires: Accurate genotype determination

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
