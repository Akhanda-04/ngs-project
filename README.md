# Bacterial Whole-Genome Sequencing (NGS) Analysis Pipeline

## Project Overview

This project implements a bacterial next-generation sequencing (NGS) analysis workflow using paired-end FASTQ reads. It covers read quality control, adapter and quality trimming, reference-based read mapping, variant calling, variant annotation, candidate prioritization, and biological interpretation.

The project uses the sequencing dataset `ERR008613` and the *Escherichia coli* K-12 substrain MG1655 reference genome (`NC_000913.3`).

The main objective is to identify genomic variants in the sequenced isolate, characterize their predicted consequences, prioritize potentially important variants, and organize the findings for further biological investigation.

> **Important:** Variant annotation and prioritization provide computational predictions. They do not independently establish a variant's biological function, pathogenicity, or antimicrobial resistance phenotype.

## Final Results Summary

The reference GFF3 originally used for annotation contained an error. It was corrected, the SnpEff database was rebuilt from the corrected annotation, and the annotation, prioritization, and interpretation stages were re-run. **All results below reflect the corrected annotation.**

| Item | Result |
|---|---|
| Annotation rows | 10 (rows, not unique variants) |
| `downstream_gene_variant` | 9 rows |
| `intergenic_region` | 1 row |
| SnpEff impact: HIGH / MODERATE / LOW | 0 / 0 / 0 |
| SnpEff impact: MODIFIER | 10 rows |
| Candidate variant | Insertion at `NC_000913.3:4296380`, `A -> ACG` |
| Candidate status | Preliminary; variant validation unresolved |

Key points:

- No HIGH-, MODERATE-, or LOW-impact consequences were predicted. Every annotation row is a MODIFIER (non-coding or neighbouring-gene) consequence.
- The only candidate is a short insertion (`A -> ACG`, i.e. a 2 bp insertion) at position 4,296,380. Its validation is unresolved.
- Because all consequences are MODIFIER-level, these results do not support any claim of altered protein function, pathogenicity, or resistance.
- Figure counts refer to annotation rows. One variant can generate several rows (for example, when it lies downstream of more than one gene), so row counts are not variant counts.

## Workflow Overview

```text
Paired-End Sequencing Data
          |
          v
   Raw FASTQ Files
          |
          v
   Initial Quality Control (FastQC, MultiQC)
          |
          v
   Read Trimming (Trimmomatic)
          |
          v
   Post-Trimming Quality Control (FastQC, MultiQC)
          |
          v
   Reference Genome Preparation
          |
          v
   Read Mapping (BWA)
          |
          v
   SAM/BAM Processing, Sorting, Indexing (SAMtools)
          |
          v
   Duplicate Marking and Alignment QC
          |
          v
   Variant Calling (GATK)
          |
          v
   GVCF and Genotyped VCF
          |
          v
   Variant Summary and QC
          |
          v
   Variant Annotation (SnpEff, corrected GFF3)
          |
          v
   Variant Prioritization
          |
          v
   IGV Visual Validation
          |
          v
   Biological Interpretation
```

## Project Directory Structure

```text
ngs/
├── data/
│   ├── raw_data/
│   └── trimmed/
├── envs/
├── ERR008613/
├── mapping/
├── qc_reports/
│   ├── pre_trim/
│   └── post_trim/
├── reference/
├── variants/
└── interpretation/
    ├── 01_variant_summary/
    ├── 02_variant_annotation/
    ├── 03_variant_prioritization/
    └── 04_biological_interpretation/
```

### Directory Description

| Directory | Purpose |
|---|---|
| `data/raw_data/` | Original paired-end FASTQ reads |
| `data/trimmed/` | Paired and unpaired reads generated after trimming |
| `envs/` | Conda environment specifications |
| `ERR008613/` | Original SRA dataset |
| `mapping/` | Alignment files, indexes, and mapping statistics |
| `qc_reports/pre_trim/` | Quality-control results before trimming |
| `qc_reports/post_trim/` | Quality-control results after trimming |
| `reference/` | Reference genome, corrected annotation files, and alignment indexes |
| `variants/` | GVCF and genotyped VCF files |
| `interpretation/01_variant_summary/` | Variant counts, genotype summaries, quality summaries, and plots |
| `interpretation/02_variant_annotation/` | SnpEff annotations and gene-level consequence summaries |
| `interpretation/03_variant_prioritization/` | Prioritized variants, candidate genes, and IGV validation records |
| `interpretation/04_biological_interpretation/` | Functional summaries, figures, and interpretation reports |

## Dataset and Reference Genome

### Sequencing dataset

- **Run accession:** `ERR008613`
- **Data type:** Paired-end sequencing reads
- **Raw read files:**
  - `data/raw_data/ERR008613_1.fastq`
  - `data/raw_data/ERR008613_2.fastq`
- **Original SRA file:** `ERR008613/ERR008613.sra`

### Reference genome

- **Organism:** *Escherichia coli*
- **Reference strain:** K-12 substrain MG1655
- **Reference accession:** `NC_000913.3`
- **Reference genome length:** Approximately 4.64 Mb
- **Reference FASTA:** `reference/reference.fasta`
- **Reference annotation (corrected):** `reference/NC_000913.3_Ecoli_MG1655.gff3`
- **GenBank annotation:** `reference/NC_000913.3_Ecoli_MG1655.gb`

The reference FASTA is indexed for alignment and downstream analysis. The reference sequence and annotation must remain consistent throughout mapping, variant calling, and annotation.

**Annotation correction:** The original GFF3 contained an error that affected annotation. The corrected GFF3 was used to rebuild the SnpEff database and regenerate all downstream annotation outputs. Any output generated before this correction should not be used.

**Reference limitation:** Using a reference genome from a different strain or lineage can introduce mapping bias and affect variant interpretation. Reference suitability should be evaluated in the context of the sample's organism and strain.

## Analysis Stages

### 1. Raw Read Quality Control

**Tools:** FastQC and MultiQC

Raw reads are assessed for sequencing quality and potential technical problems before trimming (per-base quality, per-sequence quality, GC content, adapter content, duplication, length distribution, ambiguous bases).

```text
qc_reports/pre_trim/
├── ERR008613_1_fastqc.html
├── ERR008613_2_fastqc.html
└── multiqc_report.html
```

### 2. Read Trimming

**Tool:** Trimmomatic

Adapter sequences and low-quality bases are removed. The output contains paired reads that survive trimming together and unpaired reads whose mates do not both survive.

```text
data/trimmed/
├── ERR008613_R1_paired.fastq
├── ERR008613_R2_paired.fastq
├── ERR008613_R1_unpaired.fastq
└── ERR008613_R2_unpaired.fastq
```

Compressed copies of the trimmed FASTQ files are also present.

### 3. Post-Trimming Quality Control

**Tools:** FastQC and MultiQC

```text
qc_reports/post_trim/
├── ERR008613_R1_paired_fastqc.html
├── ERR008613_R2_paired_fastqc.html
└── multiqc_report.html
```

Compare pre- and post-trimming reports to assess changes in read quality, adapter contamination, read length, and retained data.

### 4. Reference-Based Read Mapping

**Tool:** BWA

```text
reference/reference.fasta (+ .bwt, .amb, .ann, .pac, .sa)

mapping/
├── ERR008613.sam
├── ERR008613.bam
├── ERR008613.rg.bam
├── ERR008613.sorted.bam / .bai
├── ERR008613.markdup.bam / .bai
├── ERR008613_flagstats.txt
├── ERR008613_idxstats.txt
└── ERR008613.duplicate_metrics.txt
```

Mapping statistics should be reviewed to determine how many reads align and whether coverage is adequate for variant calling.

### 5. Variant Calling and Genotyping

**Tool:** GATK

```text
variants/
├── ERR008613.g.vcf.gz (+ .tbi)
├── ERR008613.haploid.vcf.gz (+ .tbi)
└── ERR008613.haploid_genotyped.vcf.gz (+ .tbi)
```

- **GVCF:** Variant-calling and reference-confidence information used in downstream genotyping.
- **Haploid VCF:** Calls produced in the haploid-calling workflow.
- **Haploid genotyped VCF:** Genotyped records used for downstream analysis.

Bacterial isolates are generally analyzed as haploid, so ploidy assumptions must match the variant-calling configuration. The exact distinction between these files should be confirmed from the commands and parameters used to generate them.

## Variant Interpretation

### Stage 1: Variant Summary and Quality Assessment

**Directory:** `interpretation/01_variant_summary/`

| File | Purpose |
|---|---|
| `variant_counts.tsv` | Overall variant counts |
| `variant_types.tsv` | Counts of variant categories |
| `quality_summary.tsv` | Summary of variant-quality metrics |
| `genotype_summary.tsv` | Sample-level genotype summary |
| `genotype_summary_per_sample.tsv` | Genotype summaries by sample |
| `titv_summary.tsv` | Transition/transversion summary |
| `run_variant_summary.sh` | Script for generating summary outputs |

Figures are in `variant_summary_figures/` (`depth_distribution.png`, `qual_distribution.png`, `quality_summary.png`, `variant_counts.png`, `variant_types.png`).

A record that lacks a filter label is not automatically equivalent to a record that passed a validated quality-filtering procedure. Interpret results in light of the filters actually applied.

### Stage 2: Variant Annotation

**Directory:** `interpretation/02_variant_annotation/`
**Tool:** SnpEff (database built from the corrected GFF3)

| File | Purpose |
|---|---|
| `ERR008613.annotated.vcf` | VCF containing SnpEff annotation |
| `annotated_variants.tsv` | Tabular variant-level annotation |
| `gene_effects.tsv` | Gene-level counts of annotated consequences |
| `snpEff_summary.csv` | Machine-readable annotation summary |
| `snpEff_summary.html` | HTML annotation report |
| `snpEff_summary.genes.txt` | Gene-oriented annotation summary |
| `snpEff_ann.log` | Annotation execution log |
| `make_annotation_tables.py` | Script for generating annotation tables |

The custom SnpEff database is stored under `interpretation/02_variant_annotation/snpeff_data/`.

**Final annotation result:** 10 annotation rows: 9 `downstream_gene_variant` and 1 `intergenic_region`, all with MODIFIER impact. No coding consequences (synonymous, missense, stop-gained, frameshift) were predicted.

Predicted consequences describe possible molecular effects. They do not by themselves prove altered gene function or a change in phenotype.

### Stage 3: Variant Prioritization and IGV Validation

**Directory:** `interpretation/03_variant_prioritization/`

| File | Purpose |
|---|---|
| `02_prioritize.sh` | Variant-prioritization workflow script |
| `high_impact_candidates.tsv` | Variants classified as high-impact by the workflow |
| `modifier_variants.tsv` | Variants assigned a MODIFIER consequence category |
| `candidate_gene_interpretation.tsv` | Candidate-gene interpretation table |
| `prioritization_notes.md` | Notes describing the prioritization approach |

Because no HIGH-impact consequences were predicted, `high_impact_candidates.tsv` is not expected to contain candidates; the MODIFIER variants are the relevant set.

#### IGV validation

```text
igv_validation/
├── igv_validation.tsv
└── screenshots/
    └── NC_000913_3_4296380.png
```

IGV allows visual inspection of reads around the candidate insertion at `NC_000913.3:4296380`. Review should assess read support, local depth, allele consistency, and possible alignment artifacts near repetitive or complex sequence.

**Status:** Validation of this candidate is unresolved. An IGV screenshot is supporting evidence for visual review, not independent proof that a variant is real or biologically meaningful.

### Stage 4: Biological Interpretation

**Directory:** `interpretation/04_biological_interpretation/`

| File | Purpose |
|---|---|
| `functional_summary.tsv` | Tabular summary of predicted functional consequences |
| `biological_findings.md` | Summary of the biological findings |
| `interpretation_report.md` | Consolidated interpretation report |
| `create_04_biological_interpretation.sh` | Script for generating interpretation outputs |
| `plot_functional_summary.py` | Script for creating functional-summary plots |

**Figures** (`interpretation/04_biological_interpretation/figures/`):

| Figure | Description |
|---|---|
| `01_annotation_consequence_counts.png` | Preliminary annotation consequence counts: 9 `downstream_gene_variant`, 1 `intergenic_region` (rows, not unique variants) |
| `02_annotation_impact_distribution.png` | Predicted impact of annotation rows: HIGH 0, MODERATE 0, LOW 0, MODIFIER 10, UNKNOWN 0 |
| `03_candidate_insertion_position.png` | Position of the candidate insertion `NC_000913.3:4296380`, `A -> ACG`, marked preliminary with validation unresolved |

The interpretation should distinguish between:

1. **Observed evidence:** Sequencing reads, called variants, and measured quality metrics.
2. **Computational predictions:** Annotated consequences and predicted molecular effects.
3. **Biological hypotheses:** Possible effects on genes, pathways, or phenotypes that require further evidence.
4. **Validated conclusions:** Claims supported by independent evidence such as experiments, reliable phenotype data, or appropriate comparative analyses.

If antimicrobial resistance is investigated, variant annotations alone are insufficient to establish resistance. Resistance claims require genotype–phenotype evidence, validated resistance determinants, and/or phenotypic susceptibility results.

## Software and Environments

```text
envs/
├── ngs-qc.yml
└── snpeff-env.from-history.yml
```

- **FastQC / MultiQC:** Read quality assessment and report aggregation
- **Trimmomatic:** Read trimming
- **BWA:** Reference-based alignment
- **SAMtools:** Alignment processing, indexing, statistics
- **GATK:** Variant calling and genotyping
- **BCFtools:** VCF inspection and processing
- **SnpEff:** Variant consequence annotation
- **IGV:** Visual inspection of candidate variants
- **Python / Bash:** Table and figure generation, workflow automation

Exact software versions and parameters should be recorded for reproducibility.

## Reproducing the Environment

```bash
conda env create -f envs/ngs-qc.yml
conda activate ngs-qc

conda env create -f envs/snpeff-env.from-history.yml
conda activate snpeff-env
```

If an environment already exists, activate it instead of creating a duplicate. Environment specifications may depend on package availability, operating-system compatibility, and the channels used when the files were created.

## Running the Analysis

Run commands from the project root:

```bash
cd "/run/media/nabil2004/New Volume/Nabil project/ngs"
```

```bash
# Generate variant-summary outputs
bash interpretation/01_variant_summary/run_variant_summary.sh

# Generate prioritized variant outputs
bash interpretation/03_variant_prioritization/02_prioritize.sh

# Generate biological interpretation outputs
bash interpretation/04_biological_interpretation/create_04_biological_interpretation.sh
```

Before running these scripts, check their input paths, environment requirements, and parameters. Re-running a script may overwrite existing outputs. If the GFF3 is changed again, rebuild the SnpEff database and re-run annotation, prioritization, and interpretation so all outputs stay consistent.

## Recommended Quality and Interpretation Checks

- [ ] Raw-read and post-trimming quality reports have been reviewed.
- [ ] Trimming parameters and retained read counts are documented.
- [ ] Mapping rate, coverage, and alignment quality have been assessed.
- [ ] The reference genome and (corrected) annotation use matching sequence identifiers.
- [ ] Variant calling and genotyping use appropriate ploidy settings.
- [ ] VCF records and variant categories were counted consistently.
- [ ] Variant-quality fields and filter status were interpreted correctly.
- [ ] SnpEff annotations were regenerated from the corrected GFF3.
- [ ] The candidate insertion at 4296380 has been reviewed in IGV and its validation resolved.
- [ ] Biological claims are separated from computational predictions.
- [ ] Software versions, scripts, and parameters are recorded.

## Limitations

- Reference bias and imperfect reference-strain matching
- Low-coverage or poorly mapped regions
- Repetitive sequences and ambiguous alignments
- False-positive and false-negative variant calls
- Limitations of automated consequence annotation
- Annotation errors caused by reference or GFF inconsistencies (one such error was found and corrected in this project)
- Very small final variant set (10 annotation rows), so conclusions are limited
- Candidate insertion not yet validated
- No independent experimental validation

The workflow does not, by itself, demonstrate causality, pathogenicity, or antimicrobial resistance.

## Future Improvements

- Resolve validation of the candidate insertion (IGV review, read support, strand balance)
- Reproducible end-to-end workflow using Snakemake or Nextflow
- Automated collection of software versions and parameters
- Coverage summaries and genome-wide coverage visualization
- More comprehensive variant filtering and sensitivity analysis
- Comparative genomics with closely related isolate genomes
- Comparison with curated antimicrobial-resistance databases when relevant
- Version control for scripts, documentation, and configuration

## Project Status

Outputs exist for quality control, read mapping, variant calling, variant summary, annotation (corrected GFF3), prioritization, IGV inspection, and biological interpretation. The existence of an output file does not by itself establish that a stage has passed validation. The current candidate variant remains **preliminary**.

## Acknowledgements

This workflow uses established open-source bioinformatics tools. Please cite the relevant software publications and databases in any formal report, thesis, or manuscript.

---

**Project:** Bacterial NGS Reference-Based Variant Analysis
**Sample/run:** `ERR008613`
**Reference:** *Escherichia coli* K-12 MG1655 (`NC_000913.3`)
**Primary outputs:** Variant summaries, annotated variants, candidate prioritization, IGV evidence, and biological interpretation reports
