#!/usr/bin/env python3
"""Merge BRAKER4 functional annotation with all-hit transcriptome TPM summary."""

import pandas as pd

annotation_file = "/work/ebg_lab/eb/diatom_consortia/functional_annotation_swissprot/06_combined_annotation/DL_diatom_master_functional_annotation.tsv"
tpm_file = "BRAKER4_summary_from_ALL_TransDecoder_ORF_hits_with_TPM.tsv"
out = "DL_diatom_BRAKER_functional_annotation_with_ALL_hit_TPM.tsv"

annot = pd.read_csv(annotation_file, sep="\t", low_memory=False)
tpm = pd.read_csv(tpm_file, sep="\t", low_memory=False)
annot.columns = annot.columns.str.strip()
tpm.columns = tpm.columns.str.strip()

possible_gene_cols = ["gene_id", "Gene_ID", "protein_id", "Protein_ID", "braker_protein_id", "query_id", "qseqid", "protein"]
gene_col = next((c for c in possible_gene_cols if c in annot.columns), annot.columns[0])
annot = annot.rename(columns={gene_col: "braker_protein_id"})

merged = annot.merge(tpm, on="braker_protein_id", how="left")
merged["compartment"] = "nuclear"
merged["average_expression_TPM"] = merged["all_hit_TPM_max"]
merged["has_transcriptome_ORF_hit"] = merged["num_unique_ORFs"].notna().map({True: "yes", False: "no"})
merged.to_csv(out, sep="\t", index=False)

print("Annotation rows:", len(annot))
print("TPM summary rows:", len(tpm))
print("Merged rows:", len(merged))
print("Genes with transcriptome ORF hit:", merged["has_transcriptome_ORF_hit"].eq("yes").sum())
print("Genes without transcriptome ORF hit:", merged["has_transcriptome_ORF_hit"].eq("no").sum())
print("Output file:", out)
