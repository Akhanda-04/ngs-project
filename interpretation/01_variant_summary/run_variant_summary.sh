
#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# NGS Variant Summary
# Input:  variants/ERR008613.haploid_genotyped.vcf.gz
# Output: interpretation/01_variant_summary/
# Requirements: bcftools, Python 3, matplotlib
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
VCF="${PROJECT_DIR}/variants/ERR008613.haploid_genotyped.vcf.gz"
OUT="${SCRIPT_DIR}"
FIG="${OUT}/variant_summary_figures"

mkdir -p "$OUT" "$FIG"

echo "============================================================"
echo "NGS Variant Summary"
echo "============================================================"
echo "Project: $PROJECT_DIR"
echo "VCF:     $VCF"
echo "Output:  $OUT"
echo

# -------------------- Check requirements ---------------------

for cmd in bcftools python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd is not installed or not in PATH." >&2
        exit 1
    fi
done

if [[ ! -s "$VCF" ]]; then
    echo "ERROR: VCF file not found or empty: $VCF" >&2
    exit 1
fi

if ! bcftools view -h "$VCF" >/dev/null; then
    echo "ERROR: Cannot read VCF: $VCF" >&2
    exit 1
fi

if ! python3 -c 'import matplotlib' >/dev/null 2>&1; then
    echo "ERROR: matplotlib is required for figures." >&2
    echo "Install with: python3 -m pip install matplotlib" >&2
    echo "Or activate an environment where matplotlib is installed." >&2
    exit 1
fi

# Validate that the VCF has at least one sample.
SAMPLES="$(bcftools query -l "$VCF")"
if [[ -z "$SAMPLES" ]]; then
    echo "ERROR: No sample columns found in VCF." >&2
    exit 1
fi

echo "Samples found:"
printf '%s\n' "$SAMPLES"
echo

# -------------------- Python analysis -----------------------

export VCF OUT FIG

python3 <<'PYTHON'
import os
import sys
import gzip
import math
import statistics
import subprocess
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VCF = os.environ["VCF"]
OUT = os.environ["OUT"]
FIG = os.environ["FIG"]

os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def write_tsv(filename, header, rows):
    path = os.path.join(OUT, filename)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")
    print(f"Created: {path}")

def fmt_number(x, digits=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        if not math.isfinite(x):
            return "NA"
        return f"{x:.{digits}f}"
    return str(x)

def parse_number(value):
    if value in (None, "", ".", "NA"):
        return None
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None

def median_or_na(values):
    return fmt_number(statistics.median(values)) if values else "NA"

def mean_or_na(values):
    return fmt_number(statistics.mean(values)) if values else "NA"

def info_dict(info_text):
    result = {}
    if info_text == ".":
        return result
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = True
    return result

def get_numeric(info, key):
    value = info.get(key)
    if value is True or value is None:
        return None
    # INFO may contain comma-separated values.
    values = str(value).split(",")
    numbers = [parse_number(v) for v in values]
    numbers = [v for v in numbers if v is not None]
    return statistics.mean(numbers) if numbers else None

def parse_format_sample(format_text, sample_text):
    keys = format_text.split(":") if format_text else []
    values = sample_text.split(":") if sample_text else []
    return dict(zip(keys, values))

def parse_gt(gt):
    if not gt or gt == ".":
        return None, None

    alleles = gt.replace("|", "/").split("/")
    if any(a == "." or a == "" for a in alleles):
        return alleles, "missing"

    try:
        allele_nums = [int(a) for a in alleles]
    except ValueError:
        return alleles, "missing"

    if all(a == 0 for a in allele_nums):
        return allele_nums, "reference"

    if len(set(allele_nums)) > 1:
        return allele_nums, "heterozygous"

    return allele_nums, "alternate"

def classify_record(ref, alts):
    """
    Mutually exclusive record-level classification.

    Multi-ALT records are counted as OTHER to avoid counting a
    single record multiple times. '*' and symbolic alleles are
    also OTHER. This is deliberately different from per-allele
    classification used by some tools.
    """
    if not alts or alts == ["."]:
        return "OTHER"

    if len(alts) != 1:
        return "OTHER"

    alt = alts[0]

    if alt == "*" or alt.startswith("<") or "[" in alt or "]" in alt:
        return "OTHER"

    if not all(c in "ACGT" for c in ref.upper() + alt.upper()):
        return "OTHER"

    ref = ref.upper()
    alt = alt.upper()

    if len(ref) == 1 and len(alt) == 1:
        return "SNP"

    if len(ref) == len(alt):
        return "MNV"

    if len(alt) > len(ref):
        return "INS"

    if len(alt) < len(ref):
        return "DEL"

    return "OTHER"

def is_transition(ref, alt):
    return (ref, alt) in {
        ("A", "G"), ("G", "A"),
        ("C", "T"), ("T", "C")
    }

# ------------------------------------------------------------
# Read VCF
# ------------------------------------------------------------

samples = []
records = 0
single_alt_records = 0
multi_alt_records = 0
star_records = 0
filter_counts = Counter()
record_types = Counter()
variant_type_records = Counter()
chrom_counts = Counter()

quals = []
info_dps = []
info_afs = []
info_mqs = []
format_dps = []
format_gqs = []
format_ads = []

ti = 0
tv = 0
titv_skipped = 0

genotype_counts = Counter()
sample_genotype_counts = defaultdict(Counter)

with gzip.open(VCF, "rt", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            samples = header[9:]
            continue

        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        chrom, pos, vid, ref, alt_text, qual_text, filt, info_text = fields[:8]
        records += 1
        chrom_counts[chrom] += 1

        alts = alt_text.split(",") if alt_text else ["."]
        if len(alts) == 1:
            single_alt_records += 1
        elif len(alts) > 1:
            multi_alt_records += 1

        if "*" in alts:
            star_records += 1

        filter_counts[filt] += 1

        # Variant class is one category per VCF record.
        vtype = classify_record(ref, alts)
        record_types[vtype] += 1
        variant_type_records[vtype] += 1

        # QUAL
        q = parse_number(qual_text)
        if q is not None:
            quals.append(q)

        # INFO fields
        info = info_dict(info_text)
        dp = get_numeric(info, "DP")
        af = get_numeric(info, "AF")
        mq = get_numeric(info, "MQ")

        if dp is not None:
            info_dps.append(dp)
        if af is not None:
            info_afs.append(af)
        if mq is not None:
            info_mqs.append(mq)

        # Ti/Tv: only unambiguous, biallelic SNPs.
        if (
            len(alts) == 1
            and len(ref) == 1
            and len(alts[0]) == 1
            and ref.upper() in "ACGT"
            and alts[0].upper() in "ACGT"
            and ref.upper() != alts[0].upper()
        ):
            if is_transition(ref.upper(), alts[0].upper()):
                ti += 1
            else:
                tv += 1
        else:
            titv_skipped += 1

        # FORMAT/sample metrics
        if len(fields) >= 9:
            fmt = fields[8]
            for sample, sample_text in zip(samples, fields[9:]):
                sample_data = parse_format_sample(fmt, sample_text)

                gt = sample_data.get("GT", ".")
                _, gt_class = parse_gt(gt)

                if gt_class is None:
                    gt_class = "missing"

                genotype_counts[gt_class] += 1
                sample_genotype_counts[sample][gt_class] += 1

                # AD: sum of allele depths when numeric.
                ad_text = sample_data.get("AD", ".")
                if ad_text not in ("", "."):
                    ad_values = [parse_number(x) for x in ad_text.split(",")]
                    ad_values = [x for x in ad_values if x is not None]
                    if ad_values:
                        format_ads.append(sum(ad_values))

                fmt_dp = parse_number(sample_data.get("DP", "."))
                fmt_gq = parse_number(sample_data.get("GQ", "."))

                if fmt_dp is not None:
                    format_dps.append(fmt_dp)
                if fmt_gq is not None:
                    format_gqs.append(fmt_gq)

# ------------------------------------------------------------
# 1. variant_counts.tsv
# ------------------------------------------------------------

count_rows = [
    ("Total variant records", records),
    ("Records with exactly one ALT allele", single_alt_records),
    ("Records with two or more ALT alleles", multi_alt_records),
    ("SNP-classified records (exclusive)", record_types["SNP"]),
    ("Insertion-classified records (exclusive)", record_types["INS"]),
    ("Deletion-classified records (exclusive)", record_types["DEL"]),
    ("MNV-classified records (exclusive)", record_types["MNV"]),
    ("OTHER-classified records (exclusive)", record_types["OTHER"]),
    ("Records containing * in ALT", star_records),
    ("Records with FILTER=PASS", filter_counts["PASS"]),
    ("Records with FILTER=.", filter_counts["."]),
    ("Records with other FILTER values",
     sum(n for k, n in filter_counts.items() if k not in ("PASS", "."))),
]
write_tsv("variant_counts.tsv", ["Metric", "Count"], count_rows)

# ------------------------------------------------------------
# 2. variant_types.tsv
# ------------------------------------------------------------

type_order = ["SNP", "INS", "DEL", "MNV", "OTHER"]
type_rows = []
for name in type_order:
    count = record_types[name]
    percentage = (100 * count / records) if records else 0
    type_rows.append((name, count, fmt_number(percentage, 2)))

write_tsv(
    "variant_types.tsv",
    ["Variant_Type", "Count", "Percentage"],
    type_rows
)

# ------------------------------------------------------------
# 3. quality_summary.tsv
# ------------------------------------------------------------

quality_metrics = [
    ("QUAL", quals),
    ("INFO_DP", info_dps),
    ("INFO_AF", info_afs),
    ("INFO_MQ", info_mqs),
    ("FORMAT_DP", format_dps),
    ("FORMAT_GQ", format_gqs),
    ("FORMAT_AD_sum", format_ads),
]

quality_rows = []
for name, values in quality_metrics:
    quality_rows.append((
        name,
        len(values),
        fmt_number(min(values)) if values else "NA",
        fmt_number(max(values)) if values else "NA",
        mean_or_na(values),
        median_or_na(values)
    ))

write_tsv(
    "quality_summary.tsv",
    ["Metric", "NonMissing_Count", "Minimum", "Maximum", "Mean", "Median"],
    quality_rows
)

# ------------------------------------------------------------
# 4. genotype_summary.tsv
# ------------------------------------------------------------

genotype_categories = [
    "reference", "alternate", "heterozygous", "missing"
]

genotype_rows = []
for category in genotype_categories:
    genotype_rows.append((
        category,
        genotype_counts[category]
    ))

genotype_rows.append((
    "Total genotype calls",
    sum(genotype_counts.values())
))

write_tsv(
    "genotype_summary.tsv",
    ["Genotype_Category", "Count"],
    genotype_rows
)

# Also report per-sample genotype counts.
per_sample_rows = []
for sample in samples:
    counts = sample_genotype_counts[sample]
    per_sample_rows.append((
        sample,
        counts["reference"],
        counts["alternate"],
        counts["heterozygous"],
        counts["missing"],
        sum(counts.values())
    ))

write_tsv(
    "genotype_summary_per_sample.tsv",
    [
        "Sample", "Reference_Calls", "Alternate_Calls",
        "Heterozygous_Calls", "Missing_Calls", "Total_Calls"
    ],
    per_sample_rows
)

# ------------------------------------------------------------
# 5. titv_summary.tsv
# ------------------------------------------------------------

total_ti_tv = ti + tv
titv_ratio = ti / tv if tv else None

titv_rows = [
    ("Transitions", ti),
    ("Transversions", tv),
    ("Ti_plus_Tv", total_ti_tv),
    ("Ti_Tv_ratio", fmt_number(titv_ratio)),
    ("Records excluded from TiTv calculation", titv_skipped),
]

write_tsv("titv_summary.tsv", ["Metric", "Value"], titv_rows)

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

def save_bar_chart(labels, values, title, ylabel, filename,
                   rotate=0):
    plt.figure(figsize=(9, 5))
    plt.bar(labels, values)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel("")
    if rotate:
        plt.xticks(rotation=rotate, ha="right")
    plt.tight_layout()
    path = os.path.join(FIG, filename)
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Created: {path}")

save_bar_chart(
    type_order,
    [record_types[x] for x in type_order],
    "Variant Types",
    "Number of VCF records",
    "variant_types.png"
)

count_plot_labels = [
    "Total records", "Single ALT", "Multi ALT",
    "SNP", "INS", "DEL", "MNV", "OTHER"
]
count_plot_values = [
    records, single_alt_records, multi_alt_records,
    record_types["SNP"], record_types["INS"],
    record_types["DEL"], record_types["MNV"],
    record_types["OTHER"]
]

save_bar_chart(
    count_plot_labels,
    count_plot_values,
    "Variant Record Summary",
    "Count",
    "variant_counts.png",
    rotate=30
)

q_labels = ["QUAL", "INFO_DP", "INFO_AF", "INFO_MQ",
            "FORMAT_DP", "FORMAT_GQ"]
q_medians = [
    statistics.median(values) if values else 0
    for _, values in quality_metrics
    if _ != "FORMAT_AD_sum"
]

save_bar_chart(
    q_labels,
    q_medians,
    "Median Quality Metrics",
    "Median value (metrics have different units)",
    "quality_summary.png",
    rotate=30
)

# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

print()
print("============================================================")
print("VARIANT SUMMARY COMPLETED")
print("============================================================")
print(f"Total VCF records: {records}")
print(f"Single-ALT records: {single_alt_records}")
print(f"Multi-ALT records: {multi_alt_records}")
print(f"SNP records: {record_types['SNP']}")
print(f"INS records: {record_types['INS']}")
print(f"DEL records: {record_types['DEL']}")
print(f"MNV records: {record_types['MNV']}")
print(f"OTHER records: {record_types['OTHER']}")
print(f"Transitions: {ti}")
print(f"Transversions: {tv}")
print(f"Ti/Tv ratio: {fmt_number(titv_ratio)}")
print(f"Samples: {len(samples)}")
print(f"Output directory: {OUT}")
print("============================================================")
PYTHON

echo
echo "Checking output files:"
find "$OUT" -maxdepth 2 -type f -printf '%P\n' | sort

echo
echo "All done."
