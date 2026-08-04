#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path.cwd()
infile = base / "07_expression/transcriptome_ORFs_vs_BRAKER4_ET_proteins.tsv"
all_hits_out = base / "07_expression/transdecoder_ORFs_vs_BRAKER4_all_hits_with_coverage.tsv"
best_out = base / "07_expression/transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv"

cols = [
    "orf_id", "braker_protein_id", "pident", "alignment_length", "orf_length_aa",
    "braker_length_aa", "orf_start", "orf_end", "braker_start", "braker_end",
    "evalue", "bitscore"
]

df = pd.read_csv(infile, sep="\t", names=cols)
for c in ["pident", "alignment_length", "orf_length_aa", "braker_length_aa", "evalue", "bitscore"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df["orf_coverage_percent"] = df["alignment_length"] / df["orf_length_aa"] * 100
df["braker_coverage_percent"] = df["alignment_length"] / df["braker_length_aa"] * 100

df = df.sort_values(
    by=["orf_id", "bitscore", "evalue", "pident", "orf_coverage_percent", "braker_coverage_percent"],
    ascending=[True, False, True, False, False, False]
)

df.to_csv(all_hits_out, sep="\t", index=False)
best = df.groupby("orf_id", as_index=False).head(1).copy()
best.to_csv(best_out, sep="\t", index=False)

print("Raw DIAMOND hit rows:", len(df))
print("Unique TransDecoder ORFs with at least one BRAKER hit:", best["orf_id"].nunique())
print("Unique BRAKER proteins hit by best ORF mappings:", best["braker_protein_id"].nunique())
print()
print("Output files:")
print(all_hits_out)
print(best_out)
