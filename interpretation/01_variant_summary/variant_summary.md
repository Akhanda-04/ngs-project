# Variant Summary

**Sample:** DRR1082025 (clinical *E. coli*, Illumina paired-end)
**Reference:** *E. coli* K-12 MG1655 (NC_000913.3)
**Caller:** GATK HaplotypeCaller → GenotypeGVCFs (default diploid ploidy)

This file summarizes the five summary tables generated from the raw, unfiltered VCF: `variant_types.tsv`, `variant_counts.tsv`, `quality_summary.tsv`, `genotype_summary.tsv`, and `genotype_categories.tsv`.

---

## 1. Variant types (`variant_types.tsv`)

Counted per ALT allele, so the total (106,945) is higher than the number of VCF records (106,801).

| Variant type | Count | Percentage |
|---|---:|---:|
| SNP | 103,256 | 96.55% |
| INS | 1,812 | 1.69% |
| DEL | 1,785 | 1.67% |
| MNV | 8 | 0.01% |
| OTHER | 84 | 0.08% |

SNPs dominate, and insertions and deletions are nearly balanced. The 84 "OTHER" alleles correspond to the 84 records with `*` (spanning deletion) in ALT.

## 2. Record-level counts (`variant_counts.tsv`)

| Metric | Count |
|---|---:|
| Total variant records | 106,801 |
| Records with exactly one ALT allele | 106,657 |
| Records with two or more ALT alleles | 144 |
| SNP-classified records | 103,229 |
| Indel-classified records | 3,590 |
| Records found in both SNP and indel sets | 18 |
| Records containing `*` in ALT | 84 |
| Records with FILTER = `.` | 106,801 |

- Multi-allelic records are rare (144, 0.13%).
- The 18 records in both the SNP and indel sets are mixed-type multi-allelic sites.
- **No filtering has been applied:** FILTER is `.` for every record.

## 3. Quality metrics (`quality_summary.tsv`)

| Metric | Min | Median | Mean | Max | n |
|---|---:|---:|---:|---:|---:|
| QUAL | 30 | 2486.06 | 2512.69 | 18809.1 | 106,801 |
| INFO_DP | 1 | 69 | 67.79 | 420 | 106,801 |
| FORMAT_DP | 1 | 67 | 66.71 | 419 | 106,801 |
| REF_AD | 0 | 0 | 0.89 | 233 | 106,801 |
| ALT_AD_total | 1 | 67 | 65.82 | 419 | 106,801 |
| AF | 0.5 | 1 | 0.991 | 1 | 106,945 |
| MQ | 21 | 60 | 59.76 | 60 | 106,801 |

- Median depth (about 67×) is consistent with the 60× mean coverage of the mapping.
- Median MQ is 60, and the median reference-supporting depth (REF_AD) is 0, as expected for homozygous-alternate calls.
- QUAL minimum of 30 is GATK's default calling threshold.
- Some calls have very low depth (minimum DP = 1) or low MQ (minimum 21), which hard filtering should remove.
- AF has 106,945 values because it is reported per ALT allele. AF = 0.5 corresponds to heterozygous genotype calls, not read-level allele fraction.

## 4. Genotype summary (`genotype_summary.tsv`)

| Metric | Count |
|---|---:|
| Total genotype records | 106,801 |
| Fully called genotypes | 106,801 |
| Missing genotypes (incl. partial) | 0 |
| Reference genotypes | 0 |
| Non-reference genotypes (fully called) | 106,801 |
| Diploid-format genotype records | 106,801 |
| Diploid homozygous alternate calls | 105,094 (98.40%) |
| Diploid heterozygous calls | 1,707 (1.60%) |
| Haploid-format / other-ploidy records | 0 |

All genotypes are fully called and all are non-reference (the VCF contains variant sites only). Every genotype is in diploid format because HaplotypeCaller was run with its default ploidy.

## 5. Genotype categories (`genotype_categories.tsv`)

| Genotype | Count | Percentage |
|---|---:|---:|
| 1/1 | 41,470 | 38.83% |
| 1\|1 | 63,624 | 59.57% |
| 0/1 | 277 | 0.26% |
| 0\|1 | 1,286 | 1.20% |
| 1/2 | 79 | 0.07% |
| 1\|2 | 65 | 0.06% |

- Homozygous alternate (`1/1` + `1|1`): 105,094 (98.40%), as expected for a haploid genome.
- Heterozygous reference/alt (`0/1` + `0|1`): 1,563 (1.46%).
- Heterozygous two-alt (`1/2` + `1|2`): 144 (0.13%), which matches the 144 multi-ALT records.
- Phased (`|`) calls are common because dense nearby variants fall on the same haplotype, which HaplotypeCaller reports as phased.

---

## Interpretation

- **Variant load:** About 103k SNPs over roughly 4.0 Mb of covered reference is about 2.6% divergence, consistent with the ~97% identity between this clinical strain and K-12. Most variants reflect lineage differences from K-12, not sequencing errors.
- **Genotype quality:** About 98.4% of calls are homozygous alternate, as expected for a haploid organism.

## Limitations

- **Ploidy:** HaplotypeCaller was run with default diploid ploidy, although bacteria are haploid. The 1,707 heterozygous calls (1.6%) are not true heterozygosity. They likely come from repetitive regions (IS elements, rRNA operons), paralogs, or misaligned reads. Future runs should use `--sample-ploidy 1`.
- **No filtering:** All 106,801 records have FILTER = `.`. Hard filtering (e.g. QD, FS, SOR, MQ, DP, alt-allele fraction) is still needed before downstream analysis.
- **Spanning deletions:** 84 records contain `*` in ALT and should be reviewed or removed after splitting multi-allelic sites.
- **Reference choice:** Reads were mapped to K-12 MG1655, so accessory-genome content of this clinical isolate (e.g. resistance genes) is not captured. Variants in highly divergent regions may be less reliable.
