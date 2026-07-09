#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

base = Path.cwd()

annotation_tpm_file = base / "07_expression/DL_diatom_master_functional_annotation_lengths_Average_TPM.tsv"

# Redo BLASTN gene-linked Phaeodactylum comparison
phaeo_hits_file = (
    base / "../../../phaeodactylum_to_diatom_blastn_redo/04_summary/"
    "phaeodactylum_vs_diatom_BLASTN_with_PT_and_BRAKER_ET_genes.tsv"
).resolve()

out_file = base / "09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv"
sorted_out_file = base / "09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_sorted_by_Average_TPM.tsv"

out_file.parent.mkdir(parents=True, exist_ok=True)


def gene_root(protein_id):
    protein_id = str(protein_id).strip()
    return re.sub(r"\.t[0-9]+$", "", protein_id)


def assign_compartment(contig):
    contig = str(contig).strip()

    plastid_contigs = {"contig_1443", "contig_4315"}
    mito_contigs = {"contig_1647", "contig_5628"}

    if contig in plastid_contigs:
        return "plastid_like"
    if contig in mito_contigs:
        return "mito_like"
    return "nuclear"


# Load annotation + TPM table
df = pd.read_csv(annotation_tpm_file, sep="\t", low_memory=False)
df.columns = df.columns.str.strip()

required_cols = [
    "protein_id",
    "transdecoder_orf_id",
    "contig_id",
    "gene_length_bp",
    "Average_TPM",
]

for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"Missing required column from annotation table: {col}")


# Assign compartment from contig ID
df["compartment"] = df["contig_id"].apply(assign_compartment)


# Load redo Phaeodactylum BLASTN gene-overlap table
phaeo = pd.read_csv(phaeo_hits_file, sep="\t", dtype=str, low_memory=False)
phaeo.columns = phaeo.columns.str.strip()

required_phaeo_cols = ["diatom_gene_id", "pt_gene_id"]

for col in required_phaeo_cols:
    if col not in phaeo.columns:
        raise ValueError(f"Missing required column from redo Phaeodactylum table: {col}")


# Keep only BLASTN hits that overlap both:
# 1) a diatom BRAKER4 ET gene model
# 2) an annotated Phaeodactylum tricornutum gene
phaeo["diatom_gene_id"] = phaeo["diatom_gene_id"].astype(str).str.strip()
phaeo["pt_gene_id"] = phaeo["pt_gene_id"].astype(str).str.strip()

valid_phaeo = phaeo[
    (~phaeo["diatom_gene_id"].isin(["", "nan", "None", "."])) &
    (~phaeo["pt_gene_id"].isin(["", "nan", "None", "."]))
].copy()

phaeo_gene_roots = set(valid_phaeo["diatom_gene_id"].dropna().astype(str))


# Match BRAKER4 isoforms to redo BLASTN gene roots
df["gene_root_for_phaeo_lookup"] = df["protein_id"].apply(gene_root)

df["in_Phaeodactylum_tricornutum"] = df["gene_root_for_phaeo_lookup"].apply(
    lambda x: "yes" if x in phaeo_gene_roots else "no"
)


# Rename BRAKER4 protein ID to final gene_id column
df = df.rename(columns={"protein_id": "gene_id"})


# Put the most important columns first
front_cols = [
    "gene_id",
    "transdecoder_orf_id",
    "compartment",
    "gene_length_bp",
    "recommended_annotation",
    "recommended_annotation_source",
    "recommended_annotation_confidence",
    "Average_TPM",
    "in_Phaeodactylum_tricornutum",
    "has_any_annotation",
    "has_swissprot_hit",
    "has_bacillariophyta_hit",
    "has_interproscan_hit",
    "has_antifam_flag",
    "antifam_accession",
    "antifam_description",
    "antifam_flag",
    "swissprot_accession",
    "swissprot_entry",
    "swissprot_protein_name",
    "swissprot_organism",
    "swissprot_gene_name",
    "swissprot_evalue",
    "swissprot_bitscore",
    "swissprot_pident",
    "swissprot_qcov_percent",
    "swissprot_scov_percent",
    "swissprot_hit_confidence",
    "bacillariophyta_accession",
    "bacillariophyta_entry",
    "bacillariophyta_protein_name",
    "bacillariophyta_organism",
    "bacillariophyta_gene_name",
    "bacillariophyta_evalue",
    "bacillariophyta_bitscore",
    "bacillariophyta_pident",
    "bacillariophyta_qcov_percent",
    "bacillariophyta_scov_percent",
    "bacillariophyta_hit_confidence",
    "interproscan_analyses",
    "signature_accessions",
    "signature_descriptions",
    "interpro_accessions",
    "interpro_descriptions",
    "go_terms",
    "pathway_annotations",
    "contig_id",
    "gene_model_start",
    "gene_model_end",
    "strand",
    "cds_length_bp",
    "protein_length_aa",
    "swissprot_taxid",
    "bacillariophyta_taxid",
    "protein_length",
    "n_interproscan_rows",
]

front_cols = [c for c in front_cols if c in df.columns]
remaining_cols = [c for c in df.columns if c not in front_cols]

final = df[front_cols + remaining_cols].copy()


# Remove temporary lookup column
final = final.drop(columns=["gene_root_for_phaeo_lookup"], errors="ignore")


# Remove detailed Phaeodactylum/BLASTN columns if they exist
pt_detail_cols = [
    "diatom_contig",
    "diatom_gene_start",
    "diatom_gene_end",
    "diatom_gene_id",
    "diatom_gene_score",
    "diatom_gene_strand",
    "hit_contig",
    "hit_start",
    "hit_end",
    "blast_hit_id",
    "bitscore",
    "blast_strand",
    "pt_contig",
    "pt_start",
    "pt_end",
    "pident",
    "aln_len",
    "evalue",
    "pt_gene_id",
    "pt_gene_name",
    "pt_gene_symbol",
    "pt_locus_tag",
    "pt_gene_strand",
]

final = final.drop(columns=[c for c in pt_detail_cols if c in final.columns], errors="ignore")


# Clean missing annotations
if "recommended_annotation" in final.columns:
    final["recommended_annotation"] = final["recommended_annotation"].fillna("unannotated")
    final.loc[
        final["recommended_annotation"].astype(str).str.strip() == "",
        "recommended_annotation"
    ] = "unannotated"


# Sort by expression
final["Average_TPM"] = pd.to_numeric(final["Average_TPM"], errors="coerce")

final.to_csv(out_file, sep="\t", index=False)

sorted_final = final.sort_values(
    by="Average_TPM",
    ascending=False,
    na_position="last"
)

sorted_final.to_csv(sorted_out_file, sep="\t", index=False)


# Summary
print("Final rows:", len(final))
print("Final columns:", len(final.columns))
print()

print("Rows with TransDecoder ORF ID:", final["transdecoder_orf_id"].notna().sum())
print("Rows without TransDecoder ORF ID:", final["transdecoder_orf_id"].isna().sum())
print()

print("Rows with Average_TPM:", final["Average_TPM"].notna().sum())
print("Rows without Average_TPM:", final["Average_TPM"].isna().sum())
print()

print("Compartment counts:")
print(final["compartment"].value_counts(dropna=False))
print()

print("Phaeodactylum yes/no counts:")
print(final["in_Phaeodactylum_tricornutum"].value_counts(dropna=False))
print()

print("Checking that no PT detail columns are present:")
bad_cols = [c for c in pt_detail_cols if c in final.columns]
print("PT detail columns found:", bad_cols)
print()

print("First 12 final columns:")
print(list(final.columns[:12]))
print()

print("Top 20 by Average_TPM:")
top_cols = [
    "gene_id",
    "transdecoder_orf_id",
    "compartment",
    "recommended_annotation",
    "Average_TPM",
    "in_Phaeodactylum_tricornutum",
]
top_cols = [c for c in top_cols if c in sorted_final.columns]
print(sorted_final[top_cols].head(20).to_string(index=False))
print()

print("Output files:")
print(out_file)
print(sorted_out_file)
