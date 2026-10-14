#!/usr/bin/env bash
# Haploid bacterial variant calling and conservative hard filtering.
# Run: THREADS=4 bash scripts/03_call_and_filter_variants.sh
# Review the filter thresholds for this dataset before interpreting results.
set -euo pipefail

SAMPLE="${SAMPLE:-DRR1082025}"
REF="${REF:-reference/reference.fasta}"
BAM="mapping/${SAMPLE}.markdup.bam"
OUT="variants"
mkdir -p "$OUT" logs

[[ -s "$REF" ]] || { echo "ERROR: Reference not found: $REF" >&2; exit 1; }
[[ -s "$BAM" ]] || { echo "ERROR: BAM not found: $BAM; run mapping script first." >&2; exit 1; }
[[ -s "${BAM}.bai" || -s "mapping/${SAMPLE}.markdup.bai" ]] || samtools index "$BAM"
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"
DICT="${REF%.*}.dict"
[[ -s "$DICT" ]] || gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

gatk HaplotypeCaller \
  -R "$REF" -I "$BAM" \
  -O "$OUT/${SAMPLE}.haploid.g.vcf.gz" \
  -ERC GVCF -ploidy 1

gatk GenotypeGVCFs \
  -R "$REF" -V "$OUT/${SAMPLE}.haploid.g.vcf.gz" \
  -O "$OUT/${SAMPLE}.haploid.vcf.gz" -ploidy 1

# Basic hard filters for a haploid bacterial callset. These are starting
# thresholds, not universal truth; inspect annotations and validate calls.
gatk VariantFiltration \
  -R "$REF" -V "$OUT/${SAMPLE}.haploid.vcf.gz" \
  --filter-name "LowQD" --filter-expression "QD < 2.0" \
  --filter-name "LowQUAL" --filter-expression "QUAL < 30.0" \
  --filter-name "LowMQ" --filter-expression "MQ < 40.0" \
  -O "$OUT/${SAMPLE}.filtered.vcf.gz"

bcftools index -f -t "$OUT/${SAMPLE}.haploid.vcf.gz"
bcftools index -f -t "$OUT/${SAMPLE}.filtered.vcf.gz"
bcftools stats "$OUT/${SAMPLE}.haploid.vcf.gz" > "$OUT/${SAMPLE}.raw.stats.txt"
bcftools stats "$OUT/${SAMPLE}.filtered.vcf.gz" > "$OUT/${SAMPLE}.filtered.stats.txt"

echo "Variant calling and filtering complete."
echo "Inspect FILTER values and validate candidate variants in IGV."
echo "Note: VariantFiltration labels failing records; it does not remove them."
echo "To create PASS-only VCF: bcftools view -f PASS -Oz -o variants/${SAMPLE}.pass.vcf.gz variants/${SAMPLE}.filtered.vcf.gz"
