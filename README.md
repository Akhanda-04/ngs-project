# NGS Read Mapping & Variant Calling

A practical bacterial NGS bioinformatics project using Illumina sequencing data to perform quality control, read trimming, reference-based read mapping, BAM processing, GATK variant calling, variant QC/filtering, and genome-browser visualization.

## Workflow

```text
SRA
 ↓
FASTQ
 ↓
FastQC → Trimmomatic → FastQC
 ↓
BWA-MEM
 ↓
SAM → BAM → sorted BAM
 ↓
Read Groups → MarkDuplicates → BAM QC
 ↓
Reference indexing (.fai + .dict)
 ↓
GATK HaplotypeCaller
 ↓
GVCF → GenotypeGVCFs → VCF
 ↓
Variant QC / Filtering
 ↓
IGV visualization
```

## Dataset

| Field | Details |
|---|---|
| SRA dataset | **DRR1082025** |
| Organism | *Escherichia coli* |
| Strain / substrain | **K-12 substr. MG1655** |
| Reference genome | **NC_000913.3** |
| Reference size | **4,641,652 bp** |
| Sequencing platform | Illumina |
| Read layout | Paired-end |
| Data source | NCBI SRA |

The reference accession **NC_000913.3** is the complete genome of *Escherichia coli* str. K-12 substr. MG1655. citeturn0search0turn0search7

> **Organism naming correction:** use **“*Escherichia coli* str. K-12 substr. MG1655”** rather than only “*E. coli* K-12” when documenting the reference genome. This makes the strain/substrain information explicit.

Large FASTQ, SAM, BAM, GVCF, and VCF data files are excluded from the repository where appropriate; small QC/index files and reports are retained.

## Results

### 1. Raw-read QC

| Metric | R1 | R2 |
|---|---:|---:|
| Reads | 1,066,446 | 1,066,446 |
| Reported sequence length | 35–301 bp | 35–301 bp |
| Per-base sequence quality | WARN | FAIL |
| Adapter Content | PASS | PASS |
| Per-sequence GC Content | FAIL | FAIL |
| Approx. GC content | ~50% | ~50% |

The raw reads showed the main QC issues documented in `pre_trim_notes.md`: poorer per-base quality in R2, localized tile-specific quality problems, GC-distribution warnings, and overrepresented technical sequences.

### 2. Mapping and BAM QC

| Metric | Result |
|---|---:|
| Total reads | 1,973,026 |
| Primary reads | 1,951,846 |
| Mapped reads | 1,474,518 |
| Mapping rate | **74.73%** |
| Primary mapped | **74.46%** |
| Properly paired | **73.95%** |
| Secondary alignments | 0 |
| Supplementary alignments | 21,180 |
| Singletons | 4,092 (0.21%) |
| Duplicate reads after MarkDuplicates | 0 |
| Mean depth | **60.09×** |
| Reference covered | **86.04%** |
| Mean mapping quality | 59 |
| Mean base quality | 36.9 |
| Zero-coverage bases | 647,951 |

The mapping QC files in the repository record the above values. The reference sequence is NC_000913.3.

### 3. Duplicate-marking metrics

Picard/GATK MarkDuplicates reported:

| Metric | Result |
|---|---:|
| Read pairs examined | 724,623 |
| Unpaired reads examined | 4,092 |
| Read-pair duplicates | 681 |
| Unpaired read duplicates | 1,020 |
| Optical duplicates | 0 |
| Percent duplication | **0.1639%** |
| Estimated library size | 385,278,606 |

Note that the `flagstat` result reports **0 duplicates** because duplicate marking was performed on the processed BAM used for that QC snapshot; the separate MarkDuplicates metrics file records the duplicate reads identified by Picard. The two files therefore should not be interpreted as identical measurements.

## Commands

The following commands document the command-line workflow used for this project. File names can be adjusted to match the local working directory.

### 1. FastQC — pre-trimming

```bash
fastqc -o qc_reports/pre_trim raw_data/DRR1082025_1.fastq.gz raw_data/DRR1082025_2.fastq.gz
```

### 2. Trimmomatic — paired-end trimming

```trimmomatic PE -threads 4 \
  raw_data/DRR1082025_1.fastq.gz \
  raw_data/DRR1082025_2.fastq.gz \
  trimmed_data/DRR1082025_1_paired.fastq.gz \
  trimmed_data/DRR1082025_1_unpaired.fastq.gz \
  trimmed_data/DRR1082025_2_paired.fastq.gz \
  trimmed_data/DRR1082025_2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
```

### 3. FastQC — post-trimming

```fastqc -o qc_reports/post_trim \
  trimmed_data/DRR1082025_1_paired.fastq.gz \
  trimmed_data/DRR1082025_2_paired.fastq.gz
```

### 4. BWA reference indexing

```bash
bwa index reference/reference.fasta
```

### 5. BWA-MEM read mapping with read groups

```bash
bwa mem -t 4 \
  -R '@RG\tID:DRR1082025\tSM:DRR1082025\tLB:Lib1\tPL:ILLUMINA\tPU:DRR1082025' \
  reference/reference.fasta \
  trimmed_data/DRR1082025_1_paired.fastq.gz \
  trimmed_data/DRR1082025_2_paired.fastq.gz \
  > mapping/DRR1082025.sam
```

### 6. SAM → BAM → sorted BAM

```bash
samtools view -@ 4 -b mapping/DRR1082025.sam > mapping/DRR1082025.bam

samtools sort -@ 4 \
  -o mapping/DRR1082025.sorted.bam \
  mapping/DRR1082025.bam

samtools index mapping/DRR1082025.sorted.bam
```

### 7. Reference indexes required by GATK

```bash
samtools faidx reference/reference.fasta

gatk CreateSequenceDictionary \
  -R reference/reference.fasta
```

This produces:

```text
reference.fasta
reference.fasta.fai
reference.dict
```

### 8. BAM QC

```bash
samtools flagstat mapping/DRR1082025.sorted.bam

samtools stats mapping/DRR1082025.sorted.bam \
  > mapping/qc/DRR1082025.stats.txt

samtools idxstats mapping/DRR1082025.sorted.bam \
  > mapping/qc/idxstats.txt

samtools coverage mapping/DRR1082025.sorted.bam \
  > mapping/qc/coverage.txt
```

### 9. Mark duplicates

```gatk MarkDuplicates \
  -I mapping/DRR1082025.rg.bam \
  -O mapping/DRR1082025.markdup.bam \
  -M mapping/DRR1082025.duplicate_metrics.txt
```

Then index the resulting BAM:

```bash
samtools index mapping/DRR1082025.markdup.bam
```

### 10. GATK HaplotypeCaller

```bash
gatk HaplotypeCaller \
  -R reference/reference.fasta \
  -I mapping/DRR1082025.markdup.bam \
  -O variants/DRR1082025.g.vcf.gz \
  -ERC GVCF
```

### 11. GenotypeGVCFs

```bash
gatk GenotypeGVCFs \
  -R reference/reference.fasta \
  -V variants/DRR1082025.g.vcf.gz \
  -O variants/DRR1082025.vcf.gz
```

### 12. IGV visualization

Load the following together in IGV:

```text
reference/reference.fasta
mapping/DRR1082025.markdup.bam
mapping/DRR1082025.markdup.bam.bai
variants/DRR1082025.vcf.gz
```

Inspect read pileups, coverage, mismatches, mapping quality, strand distribution, SNPs, and indels.

## Repository Structure

```text
ngs-project/
├── mapping/
│   ├── DRR1082025.duplicate_metrics.txt
│   └── qc/
│       ├── DRR1082025.flagstat.txt
│       ├── DRR1082025.mapping_QC.txt
│       ├── DRR1082025.stats.txt
│       ├── DRR1082025.zero_coverage_regions.txt
│       ├── coverage.txt
│       └── idxstats.txt
├── qc_reports/
│   ├── pre_trim/
│   └── post_trim/
├── reference/
│   ├── reference.fasta
│   ├── reference.fasta.fai
│   └── reference.dict
├── variants/
│   ├── DRR1082025.g.vcf.gz.tbi
│   └── DRR1082025.vcf.gz.tbi
├── pre_trim_notes.md
├── .gitignore
└── README.md
```

## Learning Goals

* NGS quality control
* Read trimming
* Reference-based mapping
* SAM/BAM processing
* Read groups
* Duplicate marking
* GATK reference preparation
* GATK HaplotypeCaller
* GVCF and VCF analysis
* Variant filtering
* IGV visualization
* Reproducible command-line bioinformatics

## Status

**Learning / Portfolio Project**

Built to develop practical skills in bacterial NGS analysis and computational genomics.

## Author

**Musaddique Mubin Nabil**

Microbiology | Bioinformatics | Computational Genomics
