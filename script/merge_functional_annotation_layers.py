#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path.cwd()
sp_file = base / "03_best_hits/DL_diatom_all_proteins_with_swissprot_annotation.tsv"
bac_file = base / "03_best_hits/DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv"
ips_file = base / "06_combined_annotation/DL_diatom_all_proteins_with_interproscan_summary.tsv"
antifam_file = base / "06_combined_annotation/DL_diatom_antifam_flagged_proteins.tsv"
outdir = base / "06_combined_annotation"
outdir.mkdir(parents=True, exist_ok=True)

out_file = outdir / "DL_diatom_master_functional_annotation.tsv"
manual_out = outdir / "DL_diatom_master_functional_annotation_for_manual_categories.tsv"
summary_out = outdir / "DL_diatom_master_functional_annotation_summary.txt"

sp = pd.read_csv(sp_file, sep="\t", dtype=str, keep_default_na=False)
bac = pd.read_csv(bac_file, sep="\t", dtype=str, keep_default_na=False)
ips = pd.read_csv(ips_file, sep="\t", dtype=str, keep_default_na=False)

for df in [sp, bac, ips]:
    df.columns = df.columns.str.strip()

merged = sp.merge(bac, on="protein_id", how="outer", suffixes=("", "_bacdup"))
merged = merged.merge(ips, on="protein_id", how="outer", suffixes=("", "_ipsdup"))

if antifam_file.exists():
    antifam = pd.read_csv(antifam_file, sep="\t", dtype=str, keep_default_na=False)
    antifam.columns = antifam.columns.str.strip()
    merged = merged.merge(antifam, on="protein_id", how="left")
else:
    merged["antifam_accession"] = ""
    merged["antifam_description"] = ""
    merged["interpretation"] = ""

merged["has_antifam_flag"] = merged.get("antifam_accession", "").astype(str).str.strip().ne("").map({True: "yes", False: "no"})
merged["antifam_flag"] = merged.get("interpretation", "")

def choose_annotation(row):
    sp_name = str(row.get("swissprot_protein_name", "")).strip()
    sp_conf = str(row.get("swissprot_hit_confidence", "")).strip()
    bac_name = str(row.get("bacillariophyta_protein_name", "")).strip()
    bac_conf = str(row.get("bacillariophyta_hit_confidence", "")).strip()
    ipr_desc = str(row.get("interpro_descriptions", "")).strip()
    sig_desc = str(row.get("signature_descriptions", "")).strip()

    if sp_name and sp_conf in {"high", "medium"}:
        return pd.Series([sp_name, "Swiss-Prot", sp_conf])
    if bac_name and bac_conf in {"high", "medium"}:
        return pd.Series([bac_name, "UniProtKB_Bacillariophyta", bac_conf])
    if sp_name:
        return pd.Series([sp_name, "Swiss-Prot", sp_conf or "reported"])
    if bac_name:
        return pd.Series([bac_name, "UniProtKB_Bacillariophyta", bac_conf or "reported"])
    if ipr_desc:
        return pd.Series([ipr_desc.split("; ")[0], "InterProScan", "domain_or_family"])
    if sig_desc:
        return pd.Series([sig_desc.split("; ")[0], "InterProScan", "signature"])
    return pd.Series(["unannotated", "none", "none"])

merged[["recommended_annotation", "recommended_annotation_source", "recommended_annotation_confidence"]] = merged.apply(choose_annotation, axis=1)

hit_cols = [c for c in ["has_swissprot_hit", "has_bacillariophyta_hit", "has_interproscan_hit"] if c in merged.columns]
for col in hit_cols:
    merged[col] = merged[col].replace({True: "yes", False: "no"}).fillna("no")
merged["has_any_annotation"] = merged[hit_cols + ["has_antifam_flag"]].eq("yes").any(axis=1).map({True: "yes", False: "no"})

front = [
    "protein_id", "recommended_annotation", "recommended_annotation_source",
    "recommended_annotation_confidence", "has_any_annotation", "has_swissprot_hit",
    "has_bacillariophyta_hit", "has_interproscan_hit", "has_antifam_flag",
    "antifam_accession", "antifam_description", "antifam_flag"
]
front = [c for c in front if c in merged.columns]
rest = [c for c in merged.columns if c not in front and not c.endswith("_bacdup") and not c.endswith("_ipsdup")]
merged = merged[front + rest]
merged.to_csv(out_file, sep="\t", index=False)
merged.to_csv(manual_out, sep="\t", index=False)

with open(summary_out, "w") as out:
    out.write(f"Rows: {len(merged)}\n")
    for col in ["has_swissprot_hit", "has_bacillariophyta_hit", "has_interproscan_hit", "has_antifam_flag", "has_any_annotation"]:
        if col in merged.columns:
            out.write(f"\n{col}\n")
            out.write(merged[col].value_counts(dropna=False).to_string())
            out.write("\n")

print("Wrote:", out_file)
print("Wrote:", manual_out)
print("Wrote:", summary_out)
