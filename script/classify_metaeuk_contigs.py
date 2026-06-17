#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path.cwd()
metaeuk_file = base / "metaeuk_output_polyp_taxonomy_tax_per_pred.tsv"
contig_to_bin_file = base / "contig_to_bin.txt"
out_csv = base / "contig_classification_final_priority.csv"

if not metaeuk_file.exists():
    raise FileNotFoundError(metaeuk_file)
if not contig_to_bin_file.exists():
    raise FileNotFoundError(contig_to_bin_file)

metaeuk = pd.read_csv(metaeuk_file, sep="\t", dtype=str, keep_default_na=False)
metaeuk.columns = metaeuk.columns.str.strip()

if "Contig_ID" not in metaeuk.columns or "Classification" not in metaeuk.columns:
    raise ValueError("Expected columns: Contig_ID and Classification")

contig_to_bin = pd.read_csv(contig_to_bin_file, sep="\t", names=["Contig_ID", "bin"], dtype=str)

def label_taxonomy(tax):
    tax = str(tax)
    if tax.strip() == "" or tax.lower() in {"nan", "none", "unclassified"}:
        return "Unclassified"
    if "d_Eukaryota" in tax or "d__Eukaryota" in tax:
        return "Eukaryota"
    if "o_Rickettsiales" in tax or "o__Rickettsiales" in tax:
        return "Mitochondria-derived"
    if "p_Cyanobacteria" in tax or "p__Cyanobacteria" in tax:
        return "Chloroplast-derived"
    if tax.strip() == "_cellular organisms":
        return "Ambiguous (Cellular Org)"
    if "d_Bacteria" in tax or "d__Bacteria" in tax:
        return "Bacteria"
    return "Other"

metaeuk["orf_label"] = metaeuk["Classification"].apply(label_taxonomy)
counts = (
    metaeuk.groupby(["Contig_ID", "orf_label"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

labels = [
    "Eukaryota",
    "Mitochondria-derived",
    "Chloroplast-derived",
    "Bacteria",
    "Ambiguous (Cellular Org)",
    "Other",
    "Unclassified",
]
for label in labels:
    if label not in counts.columns:
        counts[label] = 0

def classify(row):
    euk_group = row["Eukaryota"] + row["Mitochondria-derived"] + row["Chloroplast-derived"]
    biological = euk_group + row["Bacteria"] + row["Ambiguous (Cellular Org)"] + row["Other"]
    if biological == 0:
        return "Unclassified"
    if euk_group / biological > 0.30:
        return "Eukaryota"
    if row["Bacteria"] > euk_group and row["Bacteria"] > row["Other"]:
        return "Bacteria"
    if row["Other"] > euk_group and row["Other"] > row["Bacteria"]:
        return "Other"
    return "Ambiguous (Cellular Org)"

counts["final_contig_classification"] = counts.apply(classify, axis=1)
out = contig_to_bin.merge(counts, on="Contig_ID", how="left")
out["final_contig_classification"] = out["final_contig_classification"].fillna("Unclassified")
out.to_csv(out_csv, index=False)

print("Wrote:", out_csv)
print(out["final_contig_classification"].value_counts(dropna=False))
