#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path.cwd()
annotation_file = base / "06_combined_annotation/DL_diatom_master_functional_annotation.tsv"
gff_file = base / "01_input/DL_diatom.braker4.ET.gff3"
out_file = base / "07_expression/DL_diatom_master_functional_annotation_with_lengths.tsv"
out_file.parent.mkdir(parents=True, exist_ok=True)

def parse_attributes(attr_string):
    attrs = {}
    for item in str(attr_string).strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs

transcript_rows = []
cds_lengths = {}

with open(gff_file) as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        contig, source, feature_type, start, end, score, strand, phase, attrs_raw = parts
        start = int(start)
        end = int(end)
        attrs = parse_attributes(attrs_raw)

        if feature_type in {"mRNA", "transcript"}:
            transcript_id = attrs.get("ID", "")
            if transcript_id:
                transcript_rows.append({
                    "protein_id": transcript_id,
                    "contig_id": contig,
                    "gene_model_start": start,
                    "gene_model_end": end,
                    "strand": strand,
                    "gene_length_bp": end - start + 1,
                })
        elif feature_type == "CDS":
            parent = attrs.get("Parent", "")
            if parent:
                for transcript_id in parent.split(","):
                    cds_lengths[transcript_id] = cds_lengths.get(transcript_id, 0) + (end - start + 1)

lengths = pd.DataFrame(transcript_rows)
if lengths.empty:
    raise ValueError("No transcript or mRNA features were parsed from the GFF3.")

lengths["cds_length_bp"] = lengths["protein_id"].map(cds_lengths).fillna(0).astype(int)
lengths["protein_length_aa"] = (lengths["cds_length_bp"] // 3).astype(int)

ann = pd.read_csv(annotation_file, sep="\t", dtype=str, keep_default_na=False)
ann.columns = ann.columns.str.strip()
if "protein_id" not in ann.columns:
    raise ValueError("Expected column 'protein_id' in functional annotation table.")

merged = ann.merge(lengths, on="protein_id", how="left")
merged.to_csv(out_file, sep="\t", index=False)

print("Input annotation rows:", len(ann))
print("GFF3 transcript rows parsed:", len(lengths))
print("Output rows:", len(merged))
print()
print("Rows with gene length:", merged["gene_length_bp"].notna().sum())
print("Rows missing gene length:", merged["gene_length_bp"].isna().sum())
print()
print("Output:")
print(out_file)
