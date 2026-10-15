from Bio import SeqIO
import csv

record = SeqIO.read("reference.gbff", "genbank")

with open("gene_annotations.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow([
        "gene_name", "locus_tag", "start", "end",
        "strand", "product"
    ])

    for feature in record.features:
        if feature.type != "gene":
            continue

        q = feature.qualifiers
        writer.writerow([
            q.get("gene", [""])[0],
            q.get("locus_tag", [""])[0],
            int(feature.location.start) + 1,
            int(feature.location.end),
            feature.location.strand,
            q.get("note", [""])[0]
        ])
