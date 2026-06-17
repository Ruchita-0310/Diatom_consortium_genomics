#!/usr/bin/env python3
"""Add curation flags to the combined nuclear, plastid, and mitochondrial gene table."""

import pandas as pd

infile = "DL_diatom_final_gene_table_nuclear_plastid_mito.tsv"
outfile = "DL_diatom_final_gene_table_nuclear_plastid_mito_with_curation_flags.tsv"
curated_outfile = "DL_diatom_final_gene_table_nuclear_plastid_mito_curated_organelle_suggested.tsv"
review_outfile = "DL_diatom_organelle_rows_for_manual_review.tsv"

df = pd.read_csv(infile, sep="\t", low_memory=False)
df.columns = df.columns.str.strip()

def clean_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()

def classify_organelle_row(row):
    compartment = clean_text(row.get("compartment", ""))
    feature_type = clean_text(row.get("feature_type", ""))
    ann = clean_text(row.get("recommended_annotation", "")).lower()
    product = clean_text(row.get("product", "")).lower()
    text = " ".join([ann, product])

    if compartment == "nuclear":
        return "nuclear_row", "yes", "BRAKER nuclear prediction retained"
    if feature_type in ["tRNA", "rRNA"]:
        return "organelle_RNA", "yes", "organelle RNA feature retained"
    if feature_type == "CDS":
        if ann in ["", "unannotated"]:
            return "unannotated_CDS", "review", "unannotated organelle CDS; keep in full table, review for curated table"
        if "predicted protein" in text:
            return "low_information_CDS", "review", "generic predicted protein; review manually"
        if "hypothetical protein" in text:
            return "hypothetical_CDS", "review", "hypothetical organelle CDS; review manually"
        return "annotated_CDS", "yes", "annotated organelle CDS retained"
    return "other_feature", "review", "nonstandard feature type; review manually"

status = df.apply(classify_organelle_row, axis=1, result_type="expand")
df["curation_status"] = status[0]
df["keep_suggested"] = status[1]
df["curation_note"] = status[2]

df.to_csv(outfile, sep="\t", index=False)

curated = df[(df["compartment"] == "nuclear") | (df["compartment"].isin(["plastid", "mito"]) & df["keep_suggested"].eq("yes"))].copy()
curated.to_csv(curated_outfile, sep="\t", index=False)

review = df[df["keep_suggested"].eq("review")].copy()
review.to_csv(review_outfile, sep="\t", index=False)

print("Input rows:", len(df))
print("Full flagged output rows:", len(df))
print("Suggested curated output rows:", len(curated))
print("Review rows:", len(review))
print("Curation status counts:")
print(df["curation_status"].value_counts(dropna=False))
print("Curation status by compartment:")
print(pd.crosstab(df["compartment"], df["curation_status"]))
print("Suggested keep/review counts:")
print(pd.crosstab(df["compartment"], df["keep_suggested"]))
print("Output files:")
print(outfile)
print(curated_outfile)
print(review_outfile)
