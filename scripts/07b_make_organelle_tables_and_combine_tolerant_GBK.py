#!/usr/bin/env python3
"""
Extract plastid and mitochondrial features from GenBank files with a tolerant parser.

This parser avoids strict Biopython GenBank parsing because some GeSeq qualifiers may be malformed.
"""

import pandas as pd
import re
from collections import defaultdict

nuclear_table = "DL_diatom_final_nuclear_gene_table_for_pathway_curation.tsv"
chloroplast_gbk = "/work/ebg_lab/eb/diatom_consortia/organelle/chloro/chloroplast_contigs_core.gbk"
mito_gbk = "/work/ebg_lab/eb/diatom_consortia/organelle/mito/diatom_candidate_mitochondrion_2contigs.gbk"

plastid_out = "DL_diatom_plastid_genes_from_GBK.tsv"
mito_out = "DL_diatom_mito_genes_from_GBK.tsv"
combined_organelle_out = "DL_diatom_organelle_genes_from_GBK.tsv"
final_out = "DL_diatom_final_gene_table_nuclear_plastid_mito.tsv"

def clean_text(x):
    if x is None or pd.isna(x):
        return ""
    x = str(x).replace("\n", " ")
    return re.sub(r"\s+", " ", x).strip()

def parse_location(location):
    loc = str(location).replace(" ", "")
    strand = "-" if loc.startswith("complement") else "+"
    nums = [int(x) for x in re.findall(r"\d+", loc)]
    if not nums:
        return pd.NA, pd.NA, ".", pd.NA
    start, end = min(nums), max(nums)
    intervals = re.findall(r"(\d+)\.\.(\d+)", loc)
    if intervals:
        length = sum(abs(int(b) - int(a)) + 1 for a, b in intervals)
    elif len(nums) >= 2:
        length = abs(nums[-1] - nums[0]) + 1
    else:
        length = 1
    return start, end, strand, length

def sanitize_id(x):
    x = clean_text(x)
    x = re.sub(r"\s+", "_", x)
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    return re.sub(r"_+", "_", x).strip("_")

def parse_qualifiers(lines):
    wanted = {"gene", "product", "locus_tag", "protein_id", "note", "translation"}
    qualifiers = {}
    current_key = None
    current_value = []

    def flush():
        nonlocal current_key, current_value
        if current_key in wanted:
            val = " ".join(current_value).strip().strip('"')
            qualifiers[current_key] = clean_text(val)
        current_key = None
        current_value = []

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("/") and "=" in line:
            flush()
            key, val = line[1:].split("=", 1)
            key = key.strip()
            if key not in wanted:
                current_key = None
                current_value = []
                continue
            current_key = key
            current_value = [val.strip()]
            if current_value[0].startswith('"') and current_value[0].endswith('"') and len(current_value[0]) > 1:
                flush()
        else:
            if current_key in wanted:
                current_value.append(line)
    flush()
    return qualifiers

def read_gbk_features_tolerant(gbk_file):
    records = []
    current_contig = None
    in_features = False
    current_feature_type = None
    current_location = None
    current_qual_lines = []
    keep_types = {"CDS", "tRNA", "rRNA"}

    def flush_feature():
        nonlocal current_feature_type, current_location, current_qual_lines
        if current_feature_type in keep_types and current_location:
            records.append({
                "contig_id": current_contig,
                "feature_type": current_feature_type,
                "location": current_location,
                "qualifiers": parse_qualifiers(current_qual_lines)
            })
        current_feature_type = None
        current_location = None
        current_qual_lines = []

    with open(gbk_file, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("LOCUS"):
                flush_feature()
                parts = line.split()
                current_contig = parts[1] if len(parts) > 1 else "unknown_contig"
                in_features = False
                continue
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN") or line.startswith("//"):
                flush_feature()
                in_features = False
                continue
            if not in_features:
                continue
            m = re.match(r"^     (\S+)\s+(.+)", line)
            if m:
                flush_feature()
                current_feature_type = m.group(1)
                current_location = m.group(2).strip()
                current_qual_lines = []
                continue
            if current_feature_type:
                stripped = line.strip()
                if stripped.startswith("/"):
                    current_qual_lines.append(stripped)
                elif not current_qual_lines:
                    current_location += stripped
                else:
                    current_qual_lines.append(stripped)
    flush_feature()
    return records

def make_unique_ids(rows, compartment):
    counts = defaultdict(int)
    for row in rows:
        base = row.get("gene_id_base", "") or f"{row['contig_id']}_{row['feature_type']}_{row['feature_number']}"
        prefixed = f"{compartment}_{sanitize_id(base)}"
        counts[prefixed] += 1
        row["gene_id"] = prefixed if counts[prefixed] == 1 else f"{prefixed}_{counts[prefixed]}"
    return rows

def parse_gbk_to_table(gbk_file, compartment):
    features = read_gbk_features_tolerant(gbk_file)
    rows = []
    for i, feat in enumerate(features, start=1):
        q = feat["qualifiers"]
        start, end, strand, feature_length = parse_location(feat["location"])
        gene = clean_text(q.get("gene", ""))
        product = clean_text(q.get("product", ""))
        locus_tag = clean_text(q.get("locus_tag", ""))
        protein_id = clean_text(q.get("protein_id", ""))
        note = clean_text(q.get("note", ""))
        translation = clean_text(q.get("translation", ""))
        gene_id_base = gene or protein_id or locus_tag or sanitize_id(product.lower())[:60] or ""
        annotation = product or gene or note or "unannotated"
        if feat["feature_type"] == "CDS":
            cds_length_bp = feature_length
            protein_length_aa = len(translation.replace(" ", "")) if translation else pd.NA
        else:
            cds_length_bp = pd.NA
            protein_length_aa = pd.NA
        rows.append({
            "gene_id_base": gene_id_base, "feature_number": i, "feature_type": feat["feature_type"],
            "gene_id": "", "compartment": compartment, "gene_model_length_bp": feature_length,
            "cds_length_bp": cds_length_bp, "protein_length_aa": protein_length_aa,
            "recommended_annotation": annotation, "recommended_annotation_source": "GeSeq/GenBank",
            "recommended_annotation_confidence": "organelle_annotation", "average_expression_TPM": pd.NA,
            "all_hit_TPM_sum": pd.NA, "all_hit_TPM_max": pd.NA, "all_hit_TPM_mean": pd.NA,
            "num_ORF_hits": pd.NA, "num_unique_ORFs": pd.NA, "mapped_orf_ids": pd.NA,
            "has_transcriptome_ORF_hit": "not_tested", "in_Phaeodactylum_tricornutum": "not_tested",
            "has_any_annotation": "yes" if annotation != "unannotated" else "no",
            "has_swissprot_hit": "not_tested", "has_bacillariophyta_hit": "not_tested",
            "has_interproscan_hit": "not_tested", "has_antifam_flag": "not_tested",
            "gene_symbol": gene, "product": product, "locus_tag": locus_tag, "protein_id": protein_id,
            "note": note, "contig_id": feat["contig_id"], "gene_model_start": start,
            "gene_model_end": end, "strand": strand
        })
    return pd.DataFrame(make_unique_ids(rows, compartment)).drop(columns=["gene_id_base"], errors="ignore")

plastid = parse_gbk_to_table(chloroplast_gbk, "plastid")
mito = parse_gbk_to_table(mito_gbk, "mito")
plastid.to_csv(plastid_out, sep="\t", index=False)
mito.to_csv(mito_out, sep="\t", index=False)
organelle = pd.concat([plastid, mito], ignore_index=True)
organelle.to_csv(combined_organelle_out, sep="\t", index=False)

nuclear = pd.read_csv(nuclear_table, sep="\t", low_memory=False)
nuclear.columns = nuclear.columns.str.strip()
if "feature_type" not in nuclear.columns:
    nuclear["feature_type"] = "protein_coding_prediction"

all_cols = list(nuclear.columns)
for col in organelle.columns:
    if col not in all_cols:
        all_cols.append(col)
for col in all_cols:
    if col not in nuclear.columns:
        nuclear[col] = pd.NA
    if col not in organelle.columns:
        organelle[col] = pd.NA

combined = pd.concat([nuclear[all_cols], organelle[all_cols]], ignore_index=True)
front_cols = [
    "gene_id", "compartment", "feature_type", "gene_model_length_bp", "cds_length_bp", "protein_length_aa",
    "recommended_annotation", "recommended_annotation_source", "recommended_annotation_confidence",
    "average_expression_TPM", "in_Phaeodactylum_tricornutum", "has_transcriptome_ORF_hit",
    "has_any_annotation", "has_swissprot_hit", "has_bacillariophyta_hit", "has_interproscan_hit",
    "has_antifam_flag", "gene_symbol", "product", "locus_tag", "protein_id", "contig_id",
    "gene_model_start", "gene_model_end", "strand"
]
front_cols = [c for c in front_cols if c in combined.columns]
combined = combined[front_cols + [c for c in combined.columns if c not in front_cols]]
combined.to_csv(final_out, sep="\t", index=False)

print("Nuclear rows:", len(nuclear))
print("Plastid rows:", len(plastid))
print("Mito rows:", len(mito))
print("Combined rows:", len(combined))
print("Feature types:")
print(combined["feature_type"].value_counts(dropna=False))
print("Compartments:")
print(combined["compartment"].value_counts(dropna=False))
print("Output files:")
print(plastid_out)
print(mito_out)
print(combined_organelle_out)
print(final_out)
