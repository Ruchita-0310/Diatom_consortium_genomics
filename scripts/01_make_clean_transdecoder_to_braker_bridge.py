#!/usr/bin/env python3
"""
Parse TransDecoder ORF vs BRAKER4 DIAMOND BLASTP output.

Outputs:
  transdecoder_ORFs_vs_BRAKER4_all_hits_with_coverage.tsv
  transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv
  transdecoder_ORFs_to_BRAKER4_STRICT_for_TPM_transfer.tsv
"""

import pandas as pd

infile = "transcriptome_ORFs_vs_BRAKER4_ET_proteins.tsv"

cols = [
    "orf_id", "braker_protein_id", "pident", "alignment_length",
    "orf_length_aa", "braker_length_aa", "orf_start", "orf_end",
    "braker_start", "braker_end", "evalue", "bitscore"
]

df = pd.read_csv(infile, sep="\t", names=cols)

df["orf_coverage_percent"] = (df["alignment_length"] / df["orf_length_aa"]) * 100
df["braker_coverage_percent"] = (df["alignment_length"] / df["braker_length_aa"]) * 100
df["braker_gene_root"] = df["braker_protein_id"].str.replace(r"\.t[0-9]+$", "", regex=True)

df = df.sort_values(
    by=["orf_id", "bitscore", "evalue", "pident", "orf_coverage_percent", "braker_coverage_percent"],
    ascending=[True, False, True, False, False, False]
)

df.to_csv("transdecoder_ORFs_vs_BRAKER4_all_hits_with_coverage.tsv", sep="\t", index=False)

best = df.groupby("orf_id", as_index=False).head(1).copy()

def classify(row):
    if row["pident"] >= 95 and row["orf_coverage_percent"] >= 80:
        return "very_high_confidence"
    if row["pident"] >= 80 and row["orf_coverage_percent"] >= 70:
        return "high_confidence"
    if row["pident"] >= 50 and row["orf_coverage_percent"] >= 50:
        return "medium_confidence"
    return "low_confidence"

best["mapping_confidence"] = best.apply(classify, axis=1)
best.to_csv("transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv", sep="\t", index=False)

strict = best[(best["pident"] >= 95) & (best["orf_coverage_percent"] >= 80)].copy()
strict.to_csv("transdecoder_ORFs_to_BRAKER4_STRICT_for_TPM_transfer.tsv", sep="\t", index=False)

print("Total DIAMOND hits:", len(df))
print("ORFs with at least one BRAKER hit:", best["orf_id"].nunique())
print("BRAKER proteins hit by best ORF mappings:", best["braker_protein_id"].nunique())
print("\nMapping confidence counts:")
print(best["mapping_confidence"].value_counts())
print("\nStrict ORFs for TPM transfer:", strict["orf_id"].nunique())
print("Strict BRAKER proteins receiving TPM:", strict["braker_protein_id"].nunique())
