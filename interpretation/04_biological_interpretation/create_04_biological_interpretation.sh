#!/usr/bin/env bash
set -euo pipefail

BASE="interpretation"
ANN="$BASE/02_variant_annotation/annotated_variants.tsv"
CAND="$BASE/03_variant_prioritization/high_impact_candidates.tsv"
OUT="$BASE/04_biological_interpretation"

if [[ ! -f "$ANN" ]]; then
    echo "ERROR: Cannot find $ANN"
    exit 1
fi

mkdir -p "$OUT/figures"

python3 - "$ANN" "$CAND" "$OUT" <<'PY'
import csv
import sys
from pathlib import Path
from datetime import date

ann_path = Path(sys.argv[1])
candidate_path = Path(sys.argv[2])
out = Path(sys.argv[3])

def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t",
            lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

def val(row, *keys, default=""):
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip() not in ("", "."):
            return str(value).strip()
    return default

variants = read_tsv(ann_path)
candidates = read_tsv(candidate_path)

# Build a concise functional summary from the annotation table.
functional_fields = [
    "CHROM", "POS", "REF", "ALT", "QUAL", "DP", "GT",
    "GENE_CONTEXT", "GENE_ID", "EFFECT", "IMPACT",
    "HGVS_C", "HGVS_P", "FUNCTIONAL_INTERPRETATION",
    "EVIDENCE_STATUS"
]

functional_rows = []
for v in variants:
    effect = val(v, "effect", "EFFECT", default="unknown")
    impact = val(v, "impact", "IMPACT", default="unknown")
    gene = val(v, "gene", "GENE", default="not assigned")
    gene_id = val(v, "gene_id", "GENE_ID", default="not assigned")
    ref = val(v, "ref", "REF")
    alt = val(v, "alt", "ALT")

    if effect in ("intergenic_region", "intergenic"):
        interpretation = (
            "Located between annotated genes; direct functional effect "
            "is not established by this annotation."
        )
    elif "downstream_gene_variant" in effect:
        interpretation = (
            "Near or downstream of an annotated gene; possible regulatory "
            "or local genomic effect requires further evidence."
        )
    elif "synonymous" in effect:
        interpretation = "Synonymous coding consequence; functional effect not established."
    elif "missense" in effect:
        interpretation = "Predicted amino-acid substitution; biological effect requires evaluation."
    elif "stop_gained" in effect or "frameshift" in effect:
        interpretation = "Predicted protein-altering consequence; validate before biological conclusions."
    else:
        interpretation = (
            "Annotation describes the predicted sequence consequence; "
            "biological significance is not established."
        )

    functional_rows.append({
        "CHROM": val(v, "chrom", "CHROM"),
        "POS": val(v, "pos", "POS"),
        "REF": ref,
        "ALT": alt,
        "QUAL": val(v, "qual", "QUAL"),
        "DP": val(v, "dp", "DP"),
        "GT": val(v, "gt", "GT"),
        "GENE_CONTEXT": gene,
        "GENE_ID": gene_id,
        "EFFECT": effect,
        "IMPACT": impact,
        "HGVS_C": val(v, "hgvs_c", "HGVS_c", "HGVS_C"),
        "HGVS_P": val(v, "hgvs_p", "HGVS_p", "HGVS_P"),
        "FUNCTIONAL_INTERPRETATION": interpretation,
        "EVIDENCE_STATUS": "Preliminary; variant-calling discrepancy unresolved"
    })

write_tsv(out / "functional_summary.tsv", functional_fields, functional_rows)

# Conservative counts: annotation rows are not necessarily unique variants.
effects = [r["EFFECT"] for r in functional_rows]
impacts = [r["IMPACT"] for r in functional_rows]
n_rows = len(functional_rows)
n_high = sum(x.lower() == "high" for x in impacts)
n_moderate = sum(x.lower() == "moderate" for x in impacts)
n_intergenic = sum("intergenic" in x.lower() for x in effects)
n_downstream = sum("downstream_gene_variant" in x.lower() for x in effects)

# Findings document
finding_lines = [
    "# Biological Findings",
    "",
    f"**Generated:** {date.today().isoformat()}",
    "",
    "## Study context",
    "",
    "This analysis examines sequence variation from the bacterial sample "
    "`ERR008613` mapped against the *Escherichia coli* K-12 MG1655 reference "
    "genome (`NC_000913.3`). The findings below are based on the currently "
    "available SnpEff annotation output and are preliminary.",
    "",
    "## Current annotation findings",
    "",
    f"- Annotation-table rows: **{n_rows}**.",
    f"- Rows annotated with HIGH impact: **{n_high}**.",
    f"- Rows annotated with MODERATE impact: **{n_moderate}**.",
    f"- Rows annotated as intergenic: **{n_intergenic}**.",
    f"- Rows annotated as downstream of a gene: **{n_downstream}**.",
    "",
    "These are counts of annotation rows, not guaranteed counts of unique "
    "biological variants. One variant can have multiple transcript or gene "
    "consequences.",
    "",
    "## Interpretation",
    "",
    "The current annotation output contains a reported insertion at "
    "`NC_000913.3:4296380`, with reference allele `A` and alternate allele "
    "`ACG`. The annotation describes its context as intergenic and reports "
    "MODIFIER impact. This indicates that the current annotation does not "
    "predict a HIGH- or MODERATE-impact coding consequence for this record.",
    "",
    "The read pileup previously showed many reads with a `+2CG` insertion "
    "observation at this position. However, the GVCF-to-final-VCF discrepancy "
    "has not been resolved. Therefore, the insertion should be treated as a "
    "candidate requiring confirmation rather than a definitive biological "
    "finding.",
    "",
    "## Biological significance",
    "",
    "No conclusion about antimicrobial resistance, virulence, pathogenicity, "
    "adaptation, or a specific phenotype can be drawn from this annotation "
    "alone. An intergenic insertion could potentially affect local regulation, "
    "but that possibility requires gene-neighborhood review and independent "
    "evidence.",
    "",
    "## Limitations and next steps",
    "",
    "1. Reconcile the GVCF, genotyped VCF, and annotated VCF before final claims.",
    "2. Confirm the insertion's read-level support and allele balance.",
    "3. Review the neighboring genes and regulatory features using the reference annotation.",
    "4. Compare the candidate with validated databases or literature if a phenotype is being investigated.",
    "5. Avoid claiming resistance or functional impact without supporting evidence.",
    "",
    "## Conclusion",
    "",
    "The current output supports a preliminary annotation-level observation, "
    "not a validated biological mechanism. Findings should be updated after "
    "the variant-calling discrepancy is addressed."
]
(out / "biological_findings.md").write_text(
    "\n".join(finding_lines) + "\n", encoding="utf-8"
)

# Integrated report
report_lines = [
    "# Biological Interpretation Report",
    "",
    f"**Generated:** {date.today().isoformat()}",
    "",
    "## 1. Objective",
    "",
    "To summarize the available variant annotations and describe their "
    "potential biological relevance without overstating the evidence.",
    "",
    "## 2. Reference and sample",
    "",
    "- Sample: `ERR008613`",
    "- Reference: *Escherichia coli* K-12 MG1655",
    "- Reference contig: `NC_000913.3`",
    "- Annotation source: existing SnpEff `annotated_variants.tsv`",
    "",
    "## 3. Results",
    "",
    f"- Annotation rows available: {n_rows}",
    f"- HIGH-impact annotation rows: {n_high}",
    f"- MODERATE-impact annotation rows: {n_moderate}",
    f"- Intergenic annotation rows: {n_intergenic}",
    f"- Downstream-gene annotation rows: {n_downstream}",
    "",
    "The current annotation table reports an insertion at position "
    "`4296380` (`A` to `ACG`) as intergenic with MODIFIER impact. This is "
    "a provisional annotation and is not, by itself, evidence of a "
    "phenotypic effect.",
    "",
    "## 4. Functional interpretation",
    "",
    "Intergenic variation can sometimes influence promoters, terminators, "
    "or other regulatory elements. The present output does not establish "
    "that this insertion affects expression of a neighboring gene. No "
    "antimicrobial-resistance or virulence claim is justified from the "
    "current annotation alone.",
    "",
    "See `functional_summary.tsv` for the row-level functional summary and "
    "`biological_findings.md` for the concise interpretation.",
    "",
    "## 5. Validation status",
    "",
    "The available pileup shows strong coverage and many insertion observations "
    "at the candidate position. However, the diagnostic comparison previously "
    "performed returned no GVCF records under a filter that excluded ALT "
    "alleles containing `<NON_REF>`. That filter may exclude genuine mixed "
    "ALT records. The GVCF, genotyped VCF, and annotated VCF have therefore "
    "not yet been reconciled.",
    "",
    "The candidate is **not considered fully validated** for final biological "
    "claims. The report intentionally preserves this limitation.",
    "",
    "## 6. Limitations",
    "",
    "- One annotated record is insufficient for broad claims about genome-wide biology.",
    "- Annotation impact labels are predictions, not experimental measurements.",
    "- A single sample without a suitable comparison group cannot establish association with a phenotype.",
    "- Reference choice, variant-calling settings, ploidy, normalization, and filtering can affect the final variant set.",
    "",
    "## 7. Conclusion",
    "",
    "The current evidence supports documenting a candidate insertion and its "
    "annotated genomic context. Additional validation and biological evidence "
    "are required before assigning a functional effect or phenotype.",
    "",
    "## 8. Files",
    "",
    "- `functional_summary.tsv` — row-level functional annotation summary",
    "- `biological_findings.md` — concise biological findings",
    "- `figures/` — reserved for validated plots or evidence images",
    "",
    "## 9. Reproducibility note",
    "",
    "This report was generated from the existing annotation table. It does not "
    "rerun variant calling, independently verify the genotype, or validate "
    "the candidate against external phenotype databases."
]
(out / "interpretation_report.md").write_text(
    "\n".join(report_lines) + "\n", encoding="utf-8"
)

(out / "figures" / "README.md").write_text(
    "# Figures\n\n"
    "Place validated figures here, such as an IGV screenshot of the candidate "
    "locus or a plot generated from verified variant-summary data. Do not "
    "present a figure as validation unless its underlying evidence has been reviewed.\n",
    encoding="utf-8"
)

print(f"Created: {out}")
print(f"Annotation rows summarized: {n_rows}")
print("Files:")
for name in [
    "functional_summary.tsv",
    "biological_findings.md",
    "interpretation_report.md",
    "figures/README.md",
]:
    print(f"  {out / name}")
PY

echo
echo "========== CREATED FILES =========="
find "$OUT" -maxdepth 2 -type f -print | sort
