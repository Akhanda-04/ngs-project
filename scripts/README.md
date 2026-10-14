# Project scripts

Run all scripts from the repository root and activate the project's Conda/Mamba environment first.

## Order

1. `bash scripts/00_check_environment.sh` — check tools and create standard directories.
2. `bash scripts/01_qc_and_trim.sh` — raw FastQC, Trimmomatic, post-trim FastQC.
3. `bash scripts/02_map_and_mark_duplicates.sh` — BWA-MEM with read groups, sorted BAM, duplicate marking, and BAM QC.
4. `bash scripts/03_call_and_filter_variants.sh` — haploid GATK calling, hard-filter annotations, and bcftools summaries.

## Configuration

Scripts accept optional environment variables, for example:

```bash
THREADS=4 SAMPLE=DRR1082025 bash scripts/02_map_and_mark_duplicates.sh
```

The trimming script can use a local adapter FASTA:

```bash
ADAPTERS=/path/to/TruSeq3-PE.fa bash scripts/01_qc_and_trim.sh
```

Expected paired-read filenames are `raw_data/DRR1082025_1.fastq.gz` and `raw_data/DRR1082025_2.fastq.gz`; trimmed files are expected under `trimmed_data/`. Reference: `reference/reference.fasta`.

## Important scientific cautions

- This is a bacterial/haploid workflow; the GATK commands explicitly set ploidy to 1.
- The observed mapping rate in the current project is about 74.73%, which is unexpectedly low for the stated MG1655 reference. Investigate sample identity, read composition, and zero-coverage regions before interpreting variants.
- The filtering expressions are starting thresholds only. Review the VCF annotations and dataset characteristics; do not present filtered calls as validated without checking them.
- `gatk VariantFiltration` marks failing variants in the FILTER column; it does not remove them. Use bcftools to extract PASS records when needed.
- Large FASTQ/BAM/VCF files should remain excluded from Git. Store scripts, small metrics, and reproducible documentation instead.
- Scripts are run with `bash scripts/name.sh`; executable file mode is not required.
