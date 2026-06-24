#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

# Input files
best_hits_file = Path("blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.clean.tsv")
pt_gene_map_file = Path("blast_out/blast_hit_to_phaeodactylum_gene.tsv")

# Output file
output_file = Path("blast_out/18_diatom_nuclear_v1_best_hits_to_phaeodactylum.with_PT_genes.tsv")

# Check input files
if not best_hits_file.exists():
    raise FileNotFoundError(f"Missing input file: {best_hits_file}")

if not pt_gene_map_file.exists():
    raise FileNotFoundError(f"Missing input file: {pt_gene_map_file}")

# Read input tables
best_hits = pd.read_csv(best_hits_file, sep="\t")
pt_gene_map = pd.read_csv(pt_gene_map_file, sep="\t")

# Check required columns
required_best_cols = {"blast_hit_id"}
required_map_cols = {"blast_hit_id", "pt_gene_id", "pt_gene_name", "pt_locus_tag"}

missing_best = required_best_cols - set(best_hits.columns)
missing_map = required_map_cols - set(pt_gene_map.columns)

if missing_best:
    raise ValueError(f"Missing columns in {best_hits_file}: {missing_best}")

if missing_map:
    raise ValueError(f"Missing columns in {pt_gene_map_file}: {missing_map}")

# Some BLAST hit regions may overlap more than one Phaeodactylum gene.
# Collapse those to one row per blast_hit_id before merging.
def collapse_unique(series):
    values = (
        series.dropna()
        .astype(str)
        .replace({"": pd.NA, "NA": pd.NA})
        .dropna()
        .unique()
    )
    if len(values) == 0:
        return "NA"
    return ";".join(sorted(values))

pt_gene_map_collapsed = (
    pt_gene_map
    .groupby("blast_hit_id", as_index=False)
    .agg({
        "pt_gene_id": collapse_unique,
        "pt_gene_name": collapse_unique,
        "pt_locus_tag": collapse_unique
    })
)

# Merge Phaeodactylum gene annotations onto the clean diatom best-hit table
merged = best_hits.merge(
    pt_gene_map_collapsed,
    on="blast_hit_id",
    how="left"
)

# Fill missing Phaeodactylum gene annotations
for col in ["pt_gene_id", "pt_gene_name", "pt_locus_tag"]:
    merged[col] = merged[col].fillna("NA")

# Write output
merged.to_csv(output_file, sep="\t", index=False)

print(f"Wrote: {output_file}")
print(f"Input best-hit rows: {len(best_hits)}")
print(f"Unique BLAST hit IDs with Phaeodactylum gene overlaps: {pt_gene_map_collapsed['blast_hit_id'].nunique()}")
print(f"Output rows: {len(merged)}")
