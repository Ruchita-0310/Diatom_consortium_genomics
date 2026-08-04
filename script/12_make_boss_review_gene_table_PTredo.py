#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

base = Path.cwd()

infile = base / "09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv"

phaeo_hits_file = (
    base / "../../../phaeodactylum_to_diatom_blastn_redo/04_summary/"
    "phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv"
).resolve()

outfile = base / "09_final/DL_diatom_FINAL_gene_table_for_boss_PTredo.tsv"

def gene_root(gene_id):
    gene_id = str(gene_id).strip()
    return re.sub(r"\.t[0-9]+$", "", gene_id)

df = pd.read_csv(infile, sep="\t", low_memory=False)
df.columns = df.columns.str.strip()

required_cols = [
    "gene_id",
    "contig_id",
    "compartment",
    "gene_length_bp",
    "recommended_annotation",
    "Average_TPM",
]

for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing required column from final table: {col}")

# Rebuild Phaeodactylum yes/no from redo BLASTN table
phaeo = pd.read_csv(phaeo_hits_file, sep="\t", dtype=str, low_memory=False)
phaeo.columns = phaeo.columns.str.strip()

required_phaeo_cols = ["diatom_gene_id", "pt_gene_id"]

for col in required_phaeo_cols:
    if col not in phaeo.columns:
        raise ValueError(f"Missing required column from redo Phaeodactylum table: {col}")

phaeo["diatom_gene_id"] = phaeo["diatom_gene_id"].astype(str).str.strip()
phaeo["pt_gene_id"] = phaeo["pt_gene_id"].astype(str).str.strip()

valid_phaeo = phaeo[
    (~phaeo["diatom_gene_id"].isin(["", "nan", "None", "."])) &
    (~phaeo["pt_gene_id"].isin(["", "nan", "None", "."]))
].copy()

phaeo_gene_roots = set(valid_phaeo["diatom_gene_id"].dropna().astype(str))

df["gene_root_for_phaeo_lookup"] = df["gene_id"].apply(gene_root)

df["present_in_Phaeodactylum_tricornutum"] = df["gene_root_for_phaeo_lookup"].apply(
    lambda x: "yes" if x in phaeo_gene_roots else "no"
)

# Create simplified boss/review table
final = df[[
    "gene_id",
    "contig_id",
    "compartment",
    "gene_length_bp",
    "recommended_annotation",
    "Average_TPM",
    "present_in_Phaeodactylum_tricornutum",
]].copy()

final = final.rename(columns={
    "compartment": "diatom_compartment",
    "gene_length_bp": "diatom_gene_length_bp",
    "recommended_annotation": "functional_annotation",
    "Average_TPM": "diatom_Average_TPM",
})

final["functional_annotation"] = final["functional_annotation"].fillna("unannotated")
final.loc[
    final["functional_annotation"].astype(str).str.strip() == "",
    "functional_annotation"
] = "unannotated"

final["diatom_Average_TPM"] = pd.to_numeric(final["diatom_Average_TPM"], errors="coerce")

outfile.parent.mkdir(parents=True, exist_ok=True)
final.to_csv(outfile, sep="\t", index=False)

print("Created:", outfile)
print("Rows:", len(final))
print("Columns:", len(final.columns))
print()

print("PT redo yes/no counts:")
print(final["present_in_Phaeodactylum_tricornutum"].value_counts(dropna=False))
print()

print("Compartment counts:")
print(final["diatom_compartment"].value_counts(dropna=False))
print()

print("First 10 rows:")
print(final.head(10).to_string(index=False))
