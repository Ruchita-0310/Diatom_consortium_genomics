#!/usr/bin/env python3
"""Add gene model length, CDS length, coordinates, contig ID, and strand from BRAKER4 GFF3."""

import pandas as pd
import gzip

merged_table = "DL_diatom_BRAKER_functional_annotation_with_ALL_hit_TPM.tsv"
gff3_file = "/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/BRAKER4/final_annotation_ET/DL_diatom.braker4.ET.gff3"
out_file = "DL_diatom_BRAKER_annotation_expression_lengths.tsv"

def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs

transcript_rows = []
cds_lengths = {}

with open_maybe_gzip(gff3_file) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        seqid, source, feature, start, end, score, strand, phase, attributes = parts
        start, end = int(start), int(end)
        length = end - start + 1
        attrs = parse_attributes(attributes)
        if feature in ["mRNA", "transcript"]:
            transcript_id = attrs.get("ID", "")
            if transcript_id:
                transcript_rows.append({
                    "braker_protein_id": transcript_id,
                    "contig_id": seqid,
                    "gene_model_start": start,
                    "gene_model_end": end,
                    "strand": strand,
                    "gene_model_length_bp": length
                })
        if feature == "CDS":
            parent = attrs.get("Parent", "")
            for p in parent.split(",") if parent else []:
                cds_lengths[p] = cds_lengths.get(p, 0) + length

lengths = pd.DataFrame(transcript_rows)
cds_df = pd.DataFrame([{"braker_protein_id": k, "cds_length_bp": v} for k, v in cds_lengths.items()])
lengths = lengths.merge(cds_df, on="braker_protein_id", how="left")

df = pd.read_csv(merged_table, sep="\t", low_memory=False)
df.columns = df.columns.str.strip()
if "protein_length" in df.columns:
    df = df.rename(columns={"protein_length": "protein_length_aa"})

out = df.merge(lengths, on="braker_protein_id", how="left")
front_cols = [
    "braker_protein_id", "compartment", "gene_model_length_bp", "cds_length_bp", "protein_length_aa",
    "recommended_annotation", "recommended_annotation_source", "recommended_annotation_confidence",
    "average_expression_TPM", "all_hit_TPM_sum", "all_hit_TPM_max", "all_hit_TPM_mean",
    "num_ORF_hits", "num_unique_ORFs", "mapped_orf_ids", "has_transcriptome_ORF_hit",
    "has_any_annotation", "has_swissprot_hit", "has_bacillariophyta_hit", "has_interproscan_hit", "has_antifam_flag"
]
front_cols = [c for c in front_cols if c in out.columns]
out = out[front_cols + [c for c in out.columns if c not in front_cols]]
out.to_csv(out_file, sep="\t", index=False)

print("Input rows:", len(df))
print("Length rows parsed from GFF3:", len(lengths))
print("Output rows:", len(out))
print("Rows with gene_model_length_bp:", out["gene_model_length_bp"].notna().sum())
print("Rows missing gene_model_length_bp:", out["gene_model_length_bp"].isna().sum())
print("Rows with cds_length_bp:", out["cds_length_bp"].notna().sum())
print("Rows missing cds_length_bp:", out["cds_length_bp"].isna().sum())
print("Output file:", out_file)
