#!/usr/bin/env python3

import csv
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("interpretation/04_biological_interpretation")
INPUT = BASE / "functional_summary.tsv"
FIGURES = BASE / "figures"

FIGURES.mkdir(parents=True, exist_ok=True)

if not INPUT.exists():
    raise SystemExit(f"Input file not found: {INPUT}")

with INPUT.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

if not rows:
    raise SystemExit("The functional summary contains no rows.")

required = {"CHROM", "POS", "REF", "ALT", "EFFECT", "IMPACT"}
missing = required - set(rows[0].keys())

if missing:
    raise SystemExit(f"Missing required columns: {', '.join(sorted(missing))}")

def save_figure(fig, filename):
    path = FIGURES / filename
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Created: {path}")

# --------------------------------------------------
# FIGURE 1: Annotation consequence counts
# --------------------------------------------------
effect_counts = Counter(
    row["EFFECT"] or "Unknown" for row in rows
)

labels = [
    "\n".join(textwrap.wrap(name, width=22))
    for name in effect_counts
]
values = list(effect_counts.values())

fig, ax = plt.subplots(figsize=(9, 5.5))
bars = ax.bar(labels, values)

ax.bar_label(bars, padding=4)
ax.set_title("Preliminary Annotation Consequence Counts")
ax.set_xlabel("Predicted consequence")
ax.set_ylabel("Annotation rows (not unique variants)")
ax.set_ylim(0, max(values) * 1.2)
ax.grid(axis="y", linestyle="--", alpha=0.3)
ax.set_axisbelow(True)

fig.tight_layout()
save_figure(fig, "01_annotation_consequence_counts.png")

# --------------------------------------------------
# FIGURE 2: Annotation impact distribution
# --------------------------------------------------
impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN"]
impact_counts = Counter(
    (row["IMPACT"] or "UNKNOWN").upper() for row in rows
)

# Include standard categories even when their counts are zero.
impact_values = [impact_counts.get(name, 0) for name in impact_order]

fig, ax = plt.subplots(figsize=(8, 5.5))
bars = ax.bar(impact_order, impact_values)

ax.bar_label(bars, padding=4)
ax.set_title("Predicted Impact of Annotation Rows")
ax.set_xlabel("SnpEff impact category")
ax.set_ylabel("Annotation rows (not unique variants)")
ax.set_ylim(0, max(impact_values + [1]) * 1.25)
ax.grid(axis="y", linestyle="--", alpha=0.3)
ax.set_axisbelow(True)

fig.tight_layout()
save_figure(fig, "02_annotation_impact_distribution.png")

# --------------------------------------------------
# FIGURE 3: Candidate insertion position
# --------------------------------------------------
# All records must be checked because a single genomic
# variant can have multiple annotation rows.
variant_keys = {
    (row["CHROM"], row["POS"], row["REF"], row["ALT"])
    for row in rows
}

if len(variant_keys) != 1:
    print(
        "Skipping figure 3: the table contains "
        f"{len(variant_keys)} distinct variant representations."
    )
else:
    chrom, pos, ref, alt = next(iter(variant_keys))
    position = int(pos)

    fig, ax = plt.subplots(figsize=(10, 3.5))

    # This is a schematic position, not a read-level validation.
    ax.axhline(0, linewidth=3)
    ax.scatter([position], [0], s=110, marker="v", zorder=3)

    ax.annotate(
        f"Candidate insertion\n{chrom}:{position}\n{ref} → {alt}",
        xy=(position, 0),
        xytext=(0, 35),
        textcoords="offset points",
        ha="center",
        arrowprops={"arrowstyle": "->"},
    )

    ax.set_title("Candidate Variant Position on Reference Contig")
    ax.set_xlabel(f"Position on {chrom} (1-based bp)")
    ax.set_yticks([])
    ax.set_ylim(-0.35, 1.0)
    ax.set_xlim(position - 1000, position + 1000)
    ax.text(
        0.99, 0.03,
        "Preliminary; variant validation unresolved",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )

    fig.tight_layout()
    save_figure(fig, "03_candidate_insertion_position.png")

print("\nAnnotation consequence counts:")
for name, count in effect_counts.items():
    print(f"  {name}: {count}")

print("\nImpact counts:")
for name in impact_order:
    print(f"  {name}: {impact_counts.get(name, 0)}")

print(f"\nDistinct variant representations: {len(variant_keys)}")
print("Figure generation completed.")
