#!/usr/bin/env python3

import csv
import re
from pathlib import Path

PROJECT = Path.cwd()

VCF = PROJECT / "interpretation/02_annotation/DRR1082025.ann.vcf"
GFF = PROJECT / "reference/reference.gff3"
LOG = PROJECT / "interpretation/02_annotation/snpeff.log"

OUT = PROJECT / "interpretation/02_annotation/annotated_variants.tsv"

COLUMNS = [
    "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER",
    "DP", "AD", "AF", "REGION", "GENE_ID", "GENE_NAME",
    "GENE_STRAND", "EFFECT", "CODON_CHANGE", "PROTEIN_CHANGE",
    "ANNOTATION_SOURCE", "ANNOTATION_VERSION"
]


def parse_attributes(text):
    result = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def load_genes(gff_path):
    genes = {}

    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attr = fields

            if feature.lower() not in {"gene", "pseudogene"}:
                continue

            a = parse_attributes(attr)
            gene_name = a.get("gene", a.get("Name", ""))
            gene_id = a.get("ID", "")
            locus = a.get("locus_tag", "")

            record = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": strand,
            }

            # Index by all available stable identifiers and gene names.
            for key in (gene_id, locus, gene_name):
                if key:
                    genes.setdefault((chrom, key), record)

    return genes


def find_gene(genes, chrom, ann_gene_id, ann_gene_name):
    candidates = [
        ann_gene_id,
        ann_gene_id.removeprefix("gene-"),
        ann_gene_name,
    ]

    for candidate in candidates:
        if candidate and (chrom, candidate) in genes:
            return genes[(chrom, candidate)]

    return {}


def info_dict(info_text):
    result = {}
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def region_from_effect(effect):
    effects = set(effect.split("&"))

    if "intergenic_region" in effects:
        return "intergenic"
    if any("upstream_gene_variant" == x or
           "downstream_gene_variant" == x for x in effects):
        return "upstream/downstream"
    if "intron_variant" in effects:
        return "intronic"
    if any(x in effects for x in (
        "missense_variant", "synonymous_variant",
        "stop_gained", "stop_lost", "start_lost",
        "frameshift_variant", "inframe_insertion",
        "inframe_deletion", "protein_altering_variant",
        "splice_donor_variant", "splice_acceptor_variant"
    )):
        return "coding/splice"
    if "non_coding_transcript_exon_variant" in effects:
        return "noncoding transcript"
    if "regulatory_region_variant" in effects:
        return "regulatory"
    return "other/inspect consequence"


def first_value(value):
    if value in ("", "."):
        return ""
    return value.split(",", 1)[0]


genes = load_genes(GFF)

annotation_version = "SnpEff 5.4c; custom database ecoli_k12_custom"

rows_written = 0
records_seen = 0
records_without_ann = 0

with VCF.open() as source, OUT.open("w", newline="") as destination:
    writer = csv.DictWriter(
        destination,
        fieldnames=COLUMNS,
        delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()

    for line in source:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        records_seen += 1

        chrom, pos, record_id, ref, alt_text, qual, filt, info_text = fields[:8]
        alts = alt_text.split(",")
        info = info_dict(info_text)

        format_keys = []
        sample_values = {}

        if len(fields) >= 10:
            format_keys = fields[8].split(":")
            sample_values = dict(zip(format_keys, fields[9].split(":")))

        dp = sample_values.get("DP", info.get("DP", ""))
        ad = sample_values.get("AD", "")

        if ad == ".":
            ad = ""

        ad_values = []
        if ad:
            try:
                ad_values = [int(x) for x in ad.split(",")]
            except ValueError:
                ad_values = []

        ann_text = info.get("ANN", "")
        annotations = [x for x in ann_text.split(",") if x]

        if not annotations:
            records_without_ann += 1
            annotations = [""]

        for ann in annotations:
            ann_fields = ann.split("|") if ann else []
            ann_fields += [""] * max(0, 16 - len(ann_fields))

            allele = ann_fields[0]
            effect = ann_fields[1]
            gene_name = ann_fields[3]
            ann_gene_id = ann_fields[4]
            codon_change = ann_fields[9]
            protein_change = ann_fields[10]

            gene = find_gene(genes, chrom, ann_gene_id, gene_name)

            # Identify the alternate allele for AD-based AF calculation.
            alt_index = None
            if allele in alts:
                alt_index = alts.index(allele) + 1
            elif len(alts) == 1:
                alt_index = 1

            af = ""

            if ad_values and alt_index is not None:
                if alt_index < len(ad_values):
                    total_depth = sum(ad_values)
                    if total_depth > 0:
                        af = f"{ad_values[alt_index] / total_depth:.6f}"
            elif "AF" in sample_values and sample_values["AF"] not in ("", "."):
                af = first_value(sample_values["AF"])
            elif "AF" in info and info["AF"] not in ("", "."):
                af = first_value(info["AF"])

            writer.writerow({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt_text,
                "QUAL": qual,
                "FILTER": filt,
                "DP": dp,
                "AD": ad,
                "AF": af,
                "REGION": region_from_effect(effect) if effect else "",
                "GENE_ID": ann_gene_id or gene.get("gene_id", ""),
                "GENE_NAME": gene_name or gene.get("gene_name", ""),
                "GENE_STRAND": gene.get("strand", ""),
                "EFFECT": effect,
                "CODON_CHANGE": codon_change,
                "PROTEIN_CHANGE": protein_change,
                "ANNOTATION_SOURCE": "SnpEff",
                "ANNOTATION_VERSION": annotation_version,
            })

            rows_written += 1

print(f"Input VCF records: {records_seen}")
print(f"Output annotation rows: {rows_written}")
print(f"Records without ANN: {records_without_ann}")
print(f"Output file: {OUT}")
