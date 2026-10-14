#!/usr/bin/env bash
# Raw FastQC -> Trimmomatic paired-end trimming -> post-trim FastQC.
# Run: THREADS=4 bash scripts/01_qc_and_trim.sh
set -euo pipefail

SAMPLE="${SAMPLE:-DRR1082025}"
THREADS="${THREADS:-4}"
RAW="raw_data"
TRIM="trimmed_data"
ADAPTERS="${ADAPTERS:-}"

R1="$RAW/${SAMPLE}_1.fastq.gz"
R2="$RAW/${SAMPLE}_2.fastq.gz"
mkdir -p qc_reports/pre_trim qc_reports/post_trim "$TRIM" logs

[[ -s "$R1" && -s "$R2" ]] || {
  echo "ERROR: Expected paired FASTQs: $R1 and $R2" >&2
  echo "Download/extract the SRA data first, or set SAMPLE and paths in this script." >&2
  exit 1
}

fastqc --threads "$THREADS" -o qc_reports/pre_trim "$R1" "$R2"

# Adapter FASTA is optional; when provided, use it for ILLUMINACLIP.
CLIP=()
if [[ -n "$ADAPTERS" ]]; then
  [[ -s "$ADAPTERS" ]] || { echo "ERROR: ADAPTERS file not found: $ADAPTERS" >&2; exit 1; }
  CLIP=("ILLUMINACLIP:$ADAPTERS:2:30:10")
else
  echo "WARNING: ADAPTERS not set; trimming will use quality/length steps only."
  echo "Set ADAPTERS=/path/to/TruSeq3-PE.fa to enable adapter clipping."
fi

trimmomatic PE -threads "$THREADS" -phred33 \
  "$R1" "$R2" \
  "$TRIM/${SAMPLE}_1_paired.fastq.gz" \
  "$TRIM/${SAMPLE}_1_unpaired.fastq.gz" \
  "$TRIM/${SAMPLE}_2_paired.fastq.gz" \
  "$TRIM/${SAMPLE}_2_unpaired.fastq.gz" \
  "${CLIP[@]}" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 \
  2>&1 | tee "logs/${SAMPLE}.trimmomatic.log"

fastqc --threads "$THREADS" -o qc_reports/post_trim \
  "$TRIM/${SAMPLE}_1_paired.fastq.gz" "$TRIM/${SAMPLE}_2_paired.fastq.gz"

echo "QC and trimming complete. Review HTML reports before mapping."
