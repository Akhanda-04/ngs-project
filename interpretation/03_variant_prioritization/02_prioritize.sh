#!/usr/bin/env bash
set -euo pipefail

BASE="/run/media/nabil2004/New Volume/Nabil project/ngs"

VCF="$BASE/interpretation/02_variant_annotation/ERR008613.annotated.vcf"
OUT="$BASE/interpretation/03_variant_prioritization"

mkdir -p "$OUT"

if [[ ! -r "$VCF" ]]; then
    echo "ERROR: Cannot read VCF: $VCF" >&2
    exit 1
fi

# ---- 1. High-impact / moderate-impact candidates (per-effect rows) ----
python3 - "$VCF" "$OUT" << 'PY'
import sys, csv, collections
vcf, out = sys.argv[1], sys.argv[2]

hdr = ["CHROM","POS","REF","ALT","QUAL","DP","GT","GENE","GENE_ID","EFFECT",
       "IMPACT","HGVS_c","HGVS_p"]
rank = {"HIGH":0,"MODERATE":1,"LOW":2,"MODIFIER":3}
rows = []

with open(vcf) as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        chrom,pos,_,ref,alt,qual,flt,info = c[:8]
        gt = ""
        if len(c) > 9:
            fmt = c[8].split(":"); val = c[9].split(":")
            gt = dict(zip(fmt,val)).get("GT","")
        d = dict(kv.split("=",1) if "=" in kv else (kv,"") for kv in info.split(";"))
        
        dp = d.get("DP", "")

        if len(c) > 9:
            fmt = c[8].split(":")
            val = c[9].split(":")
            sample = dict(zip(fmt, val))
            gt = sample.get("GT", gt)
            if not dp:
                dp = sample.get("DP", "")

        if "ANN" not in d:
            continue

        annotations = []
        for ann in d["ANN"].split(","):
            a = ann.split("|")
            if len(a) >= 11:
                annotations.append(a)

        if not annotations:
            continue

        # Select the highest-impact annotation
        a = min(annotations, key=lambda x: rank.get(x[2], 9))

        rows.append([
            chrom, pos, ref, alt, qual, dp, gt,
            a[3], a[4], a[1], a[2], a[9], a[10]
        ])

rows.sort(key=lambda r:(rank.get(r[10],9), -float(r[4]) if r[4] not in (".","") else 0))

with open(f"{out}/high_impact_candidates.tsv","w",newline="") as o:
    w = csv.writer(o, delimiter="\t"); w.writerow(hdr)
    for r in rows:
        if r[10] in ("HIGH","MODERATE"): w.writerow(r)

# ---- 2. Gene-level interpretation table ----
genes = collections.defaultdict(lambda: {"id":"","n":0,"high":0,"mod":0,"effects":set(),"vars":[]})
for r in rows:
    if r[10] not in ("HIGH","MODERATE"): continue
    g = genes[r[7]]
    g["id"] = r[8]; g["n"] += 1
    g["high"] += r[10]=="HIGH"; g["mod"] += r[10]=="MODERATE"
    g["effects"].add(r[9]); g["vars"].append(f"{r[1]}:{r[2]}>{r[3]}({r[12] or r[11]})")

with open(f"{out}/candidate_gene_interpretation.tsv","w",newline="") as o:
    w = csv.writer(o, delimiter="\t")
    w.writerow(["GENE","GENE_ID","N_VARIANTS","N_HIGH","N_MODERATE","EFFECTS","VARIANTS"])
    for k,g in sorted(genes.items(), key=lambda x:(-x[1]["high"],-x[1]["n"])):
        w.writerow([k,g["id"],g["n"],g["high"],g["mod"],";".join(sorted(g["effects"])),"; ".join(g["vars"])])
PY


# ---- 2b. Save MODIFIER-impact variants for review ----
python3 - "$VCF" "$OUT" << 'PY'
import sys, csv

vcf, out = sys.argv[1], sys.argv[2]
header = [
    "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER",
    "DP", "GT", "GENE", "GENE_ID", "EFFECT",
    "IMPACT", "HGVS_c", "HGVS_p"
]
rows = []

with open(vcf) as f:
    for line in f:
        if line.startswith("#"):
            continue

        c = line.rstrip("\n").split("\t")
        if len(c) < 8:
            continue

        chrom, pos, _, ref, alt, qual, flt, info = c[:8]
        d = dict(
            kv.split("=", 1) if "=" in kv else (kv, "")
            for kv in info.split(";")
        )

        dp = d.get("DP", "")
        gt = ""

        if len(c) > 9:
            fmt = c[8].split(":")
            val = c[9].split(":")
            sample = dict(zip(fmt, val))
            gt = sample.get("GT", "")
            if not dp:
                dp = sample.get("DP", "")

        for ann in d.get("ANN", "").split(","):
            if not ann:
                continue

            a = ann.split("|")
            if len(a) < 11 or a[2] != "MODIFIER":
                continue

            rows.append([
                chrom, pos, ref, alt, qual, flt, dp, gt,
                a[3], a[4], a[1], a[2], a[9], a[10]
            ])

path = f"{out}/modifier_variants.tsv"
with open(path, "w", newline="") as o:
    w = csv.writer(o, delimiter="\t")
    w.writerow(header)
    w.writerows(rows)

print(f"MODIFIER annotations saved: {len(rows)} -> {path}")
PY

# ---- 3. Notes ----

# ---- 3. Notes ----
{
echo "# Variant prioritization notes — ERR008613"
echo
echo "- Date: $(date +%F)"
echo "- Input: $VCF"
echo "- Annotation: snpEff, E. coli K-12 MG1655, NC_000913.3"
echo "- Rule: highest-impact annotation per variant; HIGH/MODERATE candidates reported separately from MODIFIER annotations."
echo
echo "## Counts"
echo "- Total VCF records: $(grep -vc '^#' "$VCF")"
echo "- HIGH annotations: $(awk -F'\t' 'NR>1 && $11=="HIGH"' "$OUT/high_impact_candidates.tsv" | wc -l)"
echo "- MODERATE annotations: $(awk -F'\t' 'NR>1 && $11=="MODERATE"' "$OUT/high_impact_candidates.tsv" | wc -l)"
echo "- MODIFIER annotation rows: $(($(wc -l < "$OUT/modifier_variants.tsv")-1))"
echo "- Genes in HIGH/MODERATE table: $(($(wc -l < "$OUT/candidate_gene_interpretation.tsv")-1))"
echo
echo "## Effect types among HIGH/MODERATE candidates"
echo '```'
cut -f10 "$OUT/high_impact_candidates.tsv" | tail -n +2 | sort | uniq -c | sort -nr
echo '```'
echo
echo "## MODIFIER variants"
echo "See modifier_variants.tsv for downstream and intergenic annotations."
echo
echo "## Interpretation"
echo "- The current VCF contains one variant record; its annotation is MODIFIER, not HIGH or MODERATE."
echo "- The insertion at NC_000913.3:4296380 is A>ACG (insertion of CG)."
echo "- Review local read support and the nearby regulatory/intergenic region before inferring biological effects."
echo "- Confirm genotype, read depth, mapping quality, and strand-bias metrics."
} > "$OUT/prioritization_notes.md"

echo "Done. Files in $OUT:"; ls -lh "$OUT"
