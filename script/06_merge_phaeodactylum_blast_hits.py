#!/usr/bin/env python

from __future__ import print_function
import csv
from collections import OrderedDict, defaultdict

raw_blast = "02_blast/phaeodactylum_vs_diatom_dcmegablast.tsv"
pt_overlap = "03_filtered/phaeodactylum_hits_with_PT_genes.tsv"
diatom_overlap = "03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv"
out = "04_summary/phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv"

blast_cols = [
    "pt_contig", "diatom_contig", "pident", "aln_len", "mismatch", "gapopen",
    "pt_start", "pt_end", "diatom_start", "diatom_end",
    "evalue", "bitscore", "pt_len", "diatom_len", "qcovs"
]

final_cols = [
    "blast_hit_id"
] + blast_cols + [
    "pt_gene_id",
    "pt_gene_name",
    "pt_gene_symbol",
    "pt_locus_tag",
    "pt_gene_strand",
    "diatom_gene_id",
    "diatom_gene_attr_id",
    "diatom_gene_strand"
]


def clean_value(x):
    if x is None:
        return None
    x = str(x).strip()
    if x in ["", ".", "-1", "NA", "nan"]:
        return None
    return x


def add_unique(store, hit_id, field, value):
    value = clean_value(value)
    if value is None:
        return

    if hit_id not in store:
        store[hit_id] = defaultdict(list)

    if value not in store[hit_id][field]:
        store[hit_id][field].append(value)


def collapse(store, hit_id, field):
    if hit_id not in store:
        return "NA"
    vals = store[hit_id].get(field, [])
    if len(vals) == 0:
        return "NA"
    return ";".join(vals)


# ---------------------------------------------------------------------
# Read raw BLASTN output and assign stable hit IDs
# ---------------------------------------------------------------------

blast_rows = OrderedDict()

with open(raw_blast, "r") as f:
    reader = csv.reader(f, delimiter="\t")

    for i, row in enumerate(reader, start=1):
        if len(row) == 0:
            continue

        hit_id = "hit_%06d" % i
        record = OrderedDict()
        record["blast_hit_id"] = hit_id

        for col, val in zip(blast_cols, row):
            record[col] = val

        blast_rows[hit_id] = record


# ---------------------------------------------------------------------
# Read Phaeodactylum BLAST interval overlaps with PT genes
#
# Columns from:
# bedtools intersect -a phaeodactylum_blast_intervals.bed
#                    -b phaeodactylum_genes.bed
#                    -wa -wb -loj
#
# A columns: 1-13
# B columns: 14-21
# PT gene columns:
#   17 = pt_gene_id
#   18 = pt_gene_name
#   19 = pt_gene_symbol
#   20 = pt_locus_tag
#   21 = pt_gene_strand
# ---------------------------------------------------------------------

pt_store = {}

with open(pt_overlap, "r") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if len(row) < 21:
            continue

        hit_id = row[3]

        add_unique(pt_store, hit_id, "pt_gene_id", row[16])
        add_unique(pt_store, hit_id, "pt_gene_name", row[17])
        add_unique(pt_store, hit_id, "pt_gene_symbol", row[18])
        add_unique(pt_store, hit_id, "pt_locus_tag", row[19])
        add_unique(pt_store, hit_id, "pt_gene_strand", row[20])


# ---------------------------------------------------------------------
# Read diatom BLAST interval overlaps with BRAKER4 ET genes
#
# Columns from:
# bedtools intersect -a diatom_blast_intervals.bed
#                    -b diatom_BRAKER_ET_genes.bed
#                    -wa -wb -loj
#
# A columns: 1-13
# B columns: 14-19
# Diatom gene columns:
#   17 = diatom_gene_id
#   18 = diatom_gene_attr_id
#   19 = diatom_gene_strand
# ---------------------------------------------------------------------

diatom_store = {}

with open(diatom_overlap, "r") as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if len(row) < 19:
            continue

        hit_id = row[3]

        add_unique(diatom_store, hit_id, "diatom_gene_id", row[16])
        add_unique(diatom_store, hit_id, "diatom_gene_attr_id", row[17])
        add_unique(diatom_store, hit_id, "diatom_gene_strand", row[18])


# ---------------------------------------------------------------------
# Write final merged table: one row per raw BLASTN hit
# ---------------------------------------------------------------------

with open(out, "w") as f:
    writer = csv.DictWriter(f, delimiter="\t", fieldnames=final_cols, lineterminator="\n")
    writer.writeheader()

    for hit_id, record in blast_rows.items():
        record["pt_gene_id"] = collapse(pt_store, hit_id, "pt_gene_id")
        record["pt_gene_name"] = collapse(pt_store, hit_id, "pt_gene_name")
        record["pt_gene_symbol"] = collapse(pt_store, hit_id, "pt_gene_symbol")
        record["pt_locus_tag"] = collapse(pt_store, hit_id, "pt_locus_tag")
        record["pt_gene_strand"] = collapse(pt_store, hit_id, "pt_gene_strand")

        record["diatom_gene_id"] = collapse(diatom_store, hit_id, "diatom_gene_id")
        record["diatom_gene_attr_id"] = collapse(diatom_store, hit_id, "diatom_gene_attr_id")
        record["diatom_gene_strand"] = collapse(diatom_store, hit_id, "diatom_gene_strand")

        writer.writerow(record)


# ---------------------------------------------------------------------
# Print quick summary
# ---------------------------------------------------------------------

total = 0
with_pt = 0
with_diatom = 0
with_both = 0
unique_pt = set()
unique_diatom = set()

for hit_id, record in blast_rows.items():
    total += 1

    pt_genes = collapse(pt_store, hit_id, "pt_gene_id")
    diatom_genes = collapse(diatom_store, hit_id, "diatom_gene_id")

    has_pt = pt_genes != "NA"
    has_diatom = diatom_genes != "NA"

    if has_pt:
        with_pt += 1
        for g in pt_genes.split(";"):
            unique_pt.add(g)

    if has_diatom:
        with_diatom += 1
        for g in diatom_genes.split(";"):
            unique_diatom.add(g)

    if has_pt and has_diatom:
        with_both += 1

print("Wrote: %s" % out)
print("Rows: %d" % total)
print("Columns: %d" % len(final_cols))
print("total_blast_hits\t%d" % total)
print("hits_with_PT_gene\t%d" % with_pt)
print("hits_with_diatom_BRAKER_ET_gene\t%d" % with_diatom)
print("hits_with_both_PT_and_diatom_gene\t%d" % with_both)
print("unique_PT_genes_hit\t%d" % len(unique_pt))
print("unique_diatom_BRAKER_ET_genes_hit\t%d" % len(unique_diatom))
