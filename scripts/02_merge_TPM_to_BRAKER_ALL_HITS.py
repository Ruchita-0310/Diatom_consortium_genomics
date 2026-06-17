#!/usr/bin/env python3
"""
Merge transcriptome ORF average TPM values onto every TransDecoder ORF to BRAKER4 hit.

This script intentionally keeps all DIAMOND hits, not only best or strict hits.
"""

import pandas as pd

diamond_file = "transcriptome_ORFs_vs_BRAKER4_ET_proteins.tsv"
expression_file = "master_with_custom_broad_categories.csv"

diamond_cols = [
    "orf_id", "braker_protein_id", "pident", "alignment_length",
    "orf_length_aa", "braker_length_aa", "orf_start", "orf_end",
    "braker_start", "braker_end", "evalue", "bitscore"
]

hits = pd.read_csv(diamond_file, sep="\t", names=diamond_cols)
hits["orf_coverage_percent"] = (hits["alignment_length"] / hits["orf_length_aa"]) * 100
hits["braker_coverage_percent"] = (hits["alignment_length"] / hits["braker_length_aa"]) * 100
hits["braker_gene_root"] = hits["braker_protein_id"].str.replace(r"\.t[0-9]+$", "", regex=True)

expr = pd.read_csv(expression_file)
expr.columns = expr.columns.str.strip()

if "orf" not in expr.columns:
    raise ValueError("Could not find column named 'orf' in expression file.")
if "Average_TPM" not in expr.columns:
    raise ValueError("Could not find column named 'Average_TPM' in expression file.")

wanted_cols = [
    "orf", "Average_TPM", "Average_count", "TPM_sample1", "TPM_sample2", "TPM_sample3",
    "Broad_Category", "Pfam_annotation", "Pfam_Name", "ko", "ko_definition",
    "description", "Taxonomy"
]
wanted_cols = [c for c in wanted_cols if c in expr.columns]
expr_small = expr[wanted_cols].copy().rename(columns={"orf": "orf_id"})

all_hits_with_tpm = hits.merge(expr_small, on="orf_id", how="left")
all_hits_with_tpm.to_csv("ALL_DIAMOND_hits_TransDecoder_ORFs_to_BRAKER4_with_TPM.tsv", sep="\t", index=False)

all_hits_with_tpm["Average_TPM"] = pd.to_numeric(all_hits_with_tpm["Average_TPM"], errors="coerce")

def join_unique(x):
    return ";".join(sorted(set(x.dropna().astype(str))))

braker_summary = (
    all_hits_with_tpm
    .groupby("braker_protein_id", as_index=False)
    .agg(
        braker_gene_root=("braker_gene_root", "first"),
        all_hit_TPM_sum=("Average_TPM", "sum"),
        all_hit_TPM_max=("Average_TPM", "max"),
        all_hit_TPM_mean=("Average_TPM", "mean"),
        num_ORF_hits=("orf_id", "count"),
        num_unique_ORFs=("orf_id", "nunique"),
        mapped_orf_ids=("orf_id", join_unique),
        max_pident=("pident", "max"),
        max_bitscore=("bitscore", "max"),
        min_evalue=("evalue", "min"),
        max_orf_coverage_percent=("orf_coverage_percent", "max"),
        max_braker_coverage_percent=("braker_coverage_percent", "max"),
        transcriptome_descriptions=("description", join_unique),
        transcriptome_broad_categories=("Broad_Category", join_unique)
    )
)

braker_summary.to_csv("BRAKER4_summary_from_ALL_TransDecoder_ORF_hits_with_TPM.tsv", sep="\t", index=False)

print("Total DIAMOND all-hit rows:", len(hits))
print("Unique TransDecoder ORFs in all-hit file:", hits["orf_id"].nunique())
print("Unique BRAKER proteins hit:", hits["braker_protein_id"].nunique())
print("Rows after TPM merge:", len(all_hits_with_tpm))
print("Rows with Average_TPM found:", all_hits_with_tpm["Average_TPM"].notna().sum())
print("Rows missing Average_TPM:", all_hits_with_tpm["Average_TPM"].isna().sum())
print("BRAKER proteins in all-hit summary:", braker_summary["braker_protein_id"].nunique())
