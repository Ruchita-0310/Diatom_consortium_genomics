#!/usr/bin/env python3

import csv
from collections import defaultdict

clean_file = "blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.clean.tsv"
map_file = "blast_out/blast_hit_to_phaeodactylum_gene.tsv"
out_file = "blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.tsv"

hit_to_pt = defaultdict(lambda: {
    "pt_gene_id": set(),
    "pt_gene_name": set(),
    "pt_locus_tag": set()
})

with open(map_file) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        hit = row["blast_hit_id"]
        for key in ["pt_gene_id", "pt_gene_name", "pt_locus_tag"]:
            val = row[key]
            if val and val != "NA":
                hit_to_pt[hit][key].add(val)

with open(clean_file) as fin, open(out_file, "w") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    fieldnames = reader.fieldnames + ["pt_gene_id", "pt_gene_name", "pt_locus_tag"]
    writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for row in reader:
        hit = row["blast_hit_id"]
        row["pt_gene_id"] = ";".join(sorted(hit_to_pt[hit]["pt_gene_id"])) or "NA"
        row["pt_gene_name"] = ";".join(sorted(hit_to_pt[hit]["pt_gene_name"])) or "NA"
        row["pt_locus_tag"] = ";".join(sorted(hit_to_pt[hit]["pt_locus_tag"])) or "NA"
        writer.writerow(row)

print("Wrote:", out_file)
