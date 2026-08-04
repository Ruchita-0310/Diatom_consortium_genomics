#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path.cwd()
annotation_file = base / "07_expression/DL_diatom_master_functional_annotation_with_lengths.tsv"
best_mapping_file = base / "07_expression/transdecoder_ORFs_to_BRAKER4_best_hit_per_ORF.tsv"
tpm_file = base / "07_expression/master_with_custom_broad_categories.csv"
out_file = base / "07_expression/DL_diatom_master_functional_annotation_lengths_Average_TPM.tsv"

ann = pd.read_csv(annotation_file, sep="\t", dtype=str, keep_default_na=False)
best = pd.read_csv(best_mapping_file, sep="\t", low_memory=False)
tpm = pd.read_csv(tpm_file, low_memory=False)

ann.columns = ann.columns.str.strip()
best.columns = best.columns.str.strip()
tpm.columns = tpm.columns.str.strip()

required_ann = ["protein_id"]
required_best = ["orf_id", "braker_protein_id"]
required_tpm = ["orf", "Average_TPM"]
for col in required_ann:
    if col not in ann.columns:
        raise ValueError(f"Missing column in annotation table: {col}")
for col in required_best:
    if col not in best.columns:
        raise ValueError(f"Missing column in best-hit mapping table: {col}")
for col in required_tpm:
    if col not in tpm.columns:
        raise ValueError(f"Missing column in TPM file: {col}")

tpm_small = tpm[["orf", "Average_TPM"]].copy()
tpm_small["Average_TPM"] = pd.to_numeric(tpm_small["Average_TPM"], errors="coerce")

mapped = best[["orf_id", "braker_protein_id"]].merge(
    tpm_small,
    left_on="orf_id",
    right_on="orf",
    how="left",
)

# Keep only ORFs with valid numeric Average_TPM so the ORF ID and expression value stay consistent.
mapped_valid = mapped[mapped["Average_TPM"].notna()].copy()
mapped_valid = mapped_valid.sort_values(
    by=["braker_protein_id", "Average_TPM"],
    ascending=[True, False],
    na_position="last",
)

top_orf_by_braker = mapped_valid.groupby("braker_protein_id", as_index=False).head(1).copy()
top_orf_by_braker = top_orf_by_braker[["braker_protein_id", "orf_id", "Average_TPM"]].rename(
    columns={"orf_id": "transdecoder_orf_id"}
)

merged = ann.merge(
    top_orf_by_braker,
    left_on="protein_id",
    right_on="braker_protein_id",
    how="left",
)
merged = merged.drop(columns=["braker_protein_id"], errors="ignore")

old_expression_cols = [
    "average_expression_TPM", "all_hit_TPM_sum", "all_hit_TPM_max", "all_hit_TPM_mean",
    "num_ORF_hits", "num_unique_ORFs", "mapped_orf_ids", "best_hit_TPM_sum",
    "best_hit_TPM_mean", "num_best_hit_ORFs", "num_unique_best_hit_ORFs",
    "best_hit_mapped_orf_ids",
]
merged = merged.drop(columns=[c for c in old_expression_cols if c in merged.columns])

front_cols = [
    "protein_id", "transdecoder_orf_id", "recommended_annotation",
    "recommended_annotation_source", "recommended_annotation_confidence",
    "gene_length_bp", "Average_TPM",
]
front_cols = [c for c in front_cols if c in merged.columns]
remaining_cols = [c for c in merged.columns if c not in front_cols]
merged = merged[front_cols + remaining_cols]
merged.to_csv(out_file, sep="\t", index=False)

print("Annotation rows:", len(ann))
print("Best ORF-to-BRAKER rows:", len(best))
print("BRAKER proteins with Average_TPM:", merged["Average_TPM"].notna().sum())
print("BRAKER proteins without Average_TPM:", merged["Average_TPM"].isna().sum())
print("BRAKER proteins with TransDecoder ORF ID:", merged["transdecoder_orf_id"].notna().sum())
print("BRAKER proteins without TransDecoder ORF ID:", merged["transdecoder_orf_id"].isna().sum())
print("Output rows:", len(merged))
print()
print("Output:")
print(out_file)
