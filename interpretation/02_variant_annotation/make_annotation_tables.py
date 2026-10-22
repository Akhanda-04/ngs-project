#!/usr/bin/env python3

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

OUT = Path("interpretation/02_variant_annotation")
VCF = OUT / "ERR008613.annotated.vcf"
GFF = Path("reference/NC_000913.3_Ecoli_MG1655.gff3")

ANNOTATED_TSV = OUT / "annotated_variants.tsv"
GENE_TSV = OUT / "gene_effects.tsv"

ANN_FIELDS = [
    "ALLELE", "EFFECT", "IMPACT", "GENE_NAME", "GENE_ID",
    "FEATURE_TYPE", "FEATURE_ID", "TRANSCRIPT_BIOTYPE", "RANK",
    "HGVS_C", "HGVS_P", "CDNA_POSITION", "CDS_POSITION",
    "AA_POSITION", "DISTANCE", "ANNOTATION_WARNINGS"
]

# Consequences counted separately in gene_effects.tsv.
SYNONYMOUS = {"synonymous_variant"}
MISSENSE = {"missense_variant"}
STOP_GAINED = {"stop_gained"}
FRAMESHIFT = {"frameshift_variant"}

# Explicit noncoding/regulatory consequences.
NONCODING = {
    "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_region",
    "intergenic_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "regulatory_region_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "nc_transcript_variant",
    "non_coding_transcript_exon_variant",
    "gene_variant",
    "conserved_intergenic_variant",
}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_attributes(text):
    attrs = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def load_gene_metadata(gff_path):
    """Read gene coordinates and product annotations from GFF3."""
    genes = {}
    products = {}
    names_to_ids = defaultdict(list)

    with open_text(gff_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            attrs = parse_attributes(attr_text)

            feature_id = attrs.get("ID", "")
            name = (
                attrs.get("Name")
                or attrs.get("gene")
                or attrs.get("locus_tag")
                or attrs.get("gene_name")
                or ""
            )
            product = attrs.get("product", "")

            if feature.lower() == "gene" and feature_id:
                genes[feature_id] = {
                    "GENE_NAME": name,
                    "GENE_START": start,
                    "GENE_END": end,
                    "GENE_STRAND": strand,
                    "GENE_PRODUCT": product,
                }
                if name:
                    names_to_ids[name].append(feature_id)

            # Products may be stored on CDS features rather than gene features.
            if feature.lower() == "cds" and product:
                parents = attrs.get("Parent", "")
                for parent in parents.split(","):
                    if parent:
                        products[parent] = product
                if feature_id:
                    products[feature_id] = product

    for gene_id, info in genes.items():
        if not info["GENE_PRODUCT"]:
            info["GENE_PRODUCT"] = products.get(gene_id, "")

    return genes, names_to_ids


def parse_ann(info_text):
    """Return all SnpEff ANN entries from one VCF INFO column."""
    values = {}

    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            values[key] = value

    ann_text = values.get("ANN", "")
    entries = []

    for entry in ann_text.split(","):
        if not entry:
            continue
        parts = entry.split("|")
        parts += [""] * (len(ANN_FIELDS) - len(parts))
        entries.append(dict(zip(ANN_FIELDS, parts[:len(ANN_FIELDS)])))

    return entries, values


def unique_variant_key(chrom, pos, ref, alt):
    return (chrom, pos, ref, alt)


def classify_consequences(effect_text):
    effects = set(effect_text.split("&")) if effect_text else set()

    return {
        "synonymous": bool(effects & SYNONYMOUS),
        "missense": bool(effects & MISSENSE),
        "stop_gained": bool(effects & STOP_GAINED),
        "frameshift": bool(effects & FRAMESHIFT),
        "noncoding": bool(effects & NONCODING),
    }


def main():
    if not VCF.exists():
        raise SystemExit(f"ERROR: Annotated VCF not found: {VCF}")
    if not GFF.exists():
        raise SystemExit(
            f"ERROR: GFF3 not found: {GFF}\n"
            "Check the GFF3 filename under reference/ and update GFF in this script."
        )

    genes, names_to_ids = load_gene_metadata(GFF)

    # Store one set of variant keys per gene and consequence category.
    gene_variants = defaultdict(set)
    gene_counts = defaultdict(lambda: {
        "synonymous": set(),
        "missense": set(),
        "stop_gained": set(),
        "frameshift": set(),
        "noncoding": set(),
    })
    gene_names = defaultdict(set)

    annotated_columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "ANNOTATION_SOURCE", "SNPEFF_VERSION"
    ] + ANN_FIELDS

    snpeff_version = ""
    row_count = 0
    variant_count = 0

    with open_text(VCF) as source, open(ANNOTATED_TSV, "w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=annotated_columns, delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()

        for line in source:
            if line.startswith("##SnpEffVersion="):
                snpeff_version = line.strip().split("=", 1)[1].strip('"')
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, variant_id, ref, alt, qual, filt, info = fields[:8]
            ann_entries, info_values = parse_ann(info)
            variant_count += 1

            # Keep a row even if SnpEff did not report an ANN entry.
            if not ann_entries:
                ann_entries = [{}]

            for ann in ann_entries:
                row = {
                    "CHROM": chrom,
                    "POS": pos,
                    "ID": variant_id,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": filt,
                    "ANNOTATION_SOURCE": "SnpEff",
                    "SNPEFF_VERSION": snpeff_version,
                }
                row.update({field: ann.get(field, "") for field in ANN_FIELDS})
                writer.writerow(row)
                row_count += 1

                gene_id = ann.get("GENE_ID", "")
                gene_name = ann.get("GENE_NAME", "")
                effect = ann.get("EFFECT", "")
                feature_type = ann.get("FEATURE_TYPE", "")

                # Intergenic annotations describe a region between genes,
                # not an individual gene; exclude them from gene summaries.
                if (
                    not gene_id
                    or feature_type == "intergenic_region"
                    or "intergenic_region" in effect.split("&")
                    or gene_id not in genes
                ):
                    # Fall back to a gene name only if it maps uniquely.
                    candidates = names_to_ids.get(gene_name, [])
                    if len(candidates) == 1 and feature_type != "intergenic_region":
                        gene_id = candidates[0]
                    else:
                        continue

                key = unique_variant_key(chrom, pos, ref, alt)
                gene_variants[gene_id].add(key)

                if gene_name:
                    gene_names[gene_id].add(gene_name)

                categories = classify_consequences(effect)
                for category, present in categories.items():
                    if present:
                        gene_counts[gene_id][category].add(key)

    gene_columns = [
        "GENE_ID", "GENE_NAME", "GENE_PRODUCT", "GENE_START",
        "GENE_END", "GENE_STRAND", "VARIANT_COUNT",
        "SYNONYMOUS_COUNT", "MISSENSE_COUNT", "STOP_GAINED_COUNT",
        "FRAMESHIFT_COUNT", "NONCODING_VARIANT_COUNT",
        "PREDICTED_EFFECT_SUMMARY"
    ]

    with open(GENE_TSV, "w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=gene_columns, delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()

        for gene_id in sorted(gene_variants):
            info = genes[gene_id]
            counts = gene_counts[gene_id]

            name = info["GENE_NAME"] or ",".join(sorted(gene_names[gene_id]))
            summary = []

            for label, category in [
                ("synonymous", "synonymous"),
                ("missense", "missense"),
                ("stop_gained", "stop_gained"),
                ("frameshift", "frameshift"),
                ("noncoding", "noncoding"),
            ]:
                count = len(counts[category])
                if count:
                    summary.append(f"{label}={count}")

            writer.writerow({
                "GENE_ID": gene_id,
                "GENE_NAME": name,
                "GENE_PRODUCT": info["GENE_PRODUCT"],
                "GENE_START": info["GENE_START"],
                "GENE_END": info["GENE_END"],
                "GENE_STRAND": info["GENE_STRAND"],
                "VARIANT_COUNT": len(gene_variants[gene_id]),
                "SYNONYMOUS_COUNT": len(counts["synonymous"]),
                "MISSENSE_COUNT": len(counts["missense"]),
                "STOP_GAINED_COUNT": len(counts["stop_gained"]),
                "FRAMESHIFT_COUNT": len(counts["frameshift"]),
                "NONCODING_VARIANT_COUNT": len(counts["noncoding"]),
                "PREDICTED_EFFECT_SUMMARY": ";".join(summary),
            })

    print(f"Annotated VCF records processed: {variant_count}")
    print(f"Individual annotation rows written: {row_count}")
    print(f"Genes summarized: {len(gene_variants)}")
    print(f"Created: {ANNOTATED_TSV}")
    print(f"Created: {GENE_TSV}")
    print("Note: consequence counts are unique CHROM/POS/REF/ALT variants per gene.")


if __name__ == "__main__":
    main()
