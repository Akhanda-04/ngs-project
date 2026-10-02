# NGS Read Mapping & Variant Calling

A practical bacterial NGS bioinformatics project using Illumina sequencing data to perform quality control, read mapping, BAM processing, variant calling, filtering, and visualization.

## Workflow

```text
SRA
 ↓
FASTQ
 ↓
FastQC → Trimmomatic
 ↓
BWA
 ↓
SAM → BAM
 ↓
SAMtools
 ↓
GATK HaplotypeCaller
 ↓
GVCF → VCF
 ↓
Variant Filtering
 ↓
IGV
```

## Dataset

* **Dataset:** DRR1082025
* **Organism:** *Escherichia coli* K-12
* **Reference:** NC_000913.3
* **Platform:** Illumina
* **Source:** NCBI SRA

Large sequencing and alignment files are excluded from this repository.

## Tools

* FastQC
* Trimmomatic
* BWA
* SAMtools
* GATK
* IGV
* SRA Toolkit
* Conda/Mamba
* Linux

## Repository Structure

```text
ngs-project/
├── mapping/
├── qc_reports/
├── reference/
├── variants/
├── pre_trim_notes.md
├── .gitignore
└── README.md
```

## Learning Goals

* NGS quality control
* Read trimming
* Reference-based mapping
* SAM/BAM processing
* GATK variant calling
* VCF analysis
* Variant filtering
* IGV visualization

## Status

**Learning / Portfolio Project**

Built to develop practical skills in bacterial NGS analysis and computational genomics.

## Author

**Akhanda-04**

Microbiology | Bioinformatics | Computational Genomics
