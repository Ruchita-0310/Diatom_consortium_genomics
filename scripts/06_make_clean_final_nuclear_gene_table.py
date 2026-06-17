#!/usr/bin/env python3
"""Create a clean nuclear gene table for pathway curation."""

import pandas as pd

infile = "DL_diatom_BRAKER_annotation_expression_lengths_ALL_Phaeodactylum_hits.tsv"
outfile = "DL_diatom_final_nuclear_gene_table_for_pathway_curation.tsv"

df = pd.read_csv(infile, sep="\t", low_memory=False)
df.columns = df.columns.str.strip()

keep_cols = [
    "braker_protein_id", "compartment", "gene_model_length_bp", "cds_length_bp", "protein_length_aa",
    "recommended_annotation", "recommended_annotation_source", "recommended_annotation_confidence",
    "average_expression_TPM", "all_hit_TPM_sum", "all_hit_TPM_max", "all_hit_TPM_mean",
    "num_ORF_hits", "num_unique_ORFs", "mapped_orf_ids", "has_transcriptome_ORF_hit",
    "in_Phaeodactylum_tricornutum", "phaeodactylum_num_all_hits", "phaeodactylum_num_unique_PT_genes",
    "phaeodactylum_all_pt_gene_ids", "phaeodactylum_all_pt_gene_names", "phaeodactylum_all_pt_locus_tags",
    "phaeodactylum_best_bitscore", "phaeodactylum_best_pident", "phaeodactylum_best_evalue",
    "phaeodactylum_max_alignment_length", "has_any_annotation", "has_swissprot_hit", "has_bacillariophyta_hit",
    "has_interproscan_hit", "has_antifam_flag", "swissprot_accession", "swissprot_entry",
    "swissprot_protein_name", "swissprot_organism", "swissprot_gene_name", "swissprot_evalue",
    "swissprot_bitscore", "swissprot_pident", "swissprot_qcov_percent", "swissprot_scov_percent",
    "swissprot_hit_confidence", "bacillariophyta_accession", "bacillariophyta_entry",
    "bacillariophyta_protein_name", "bacillariophyta_organism", "bacillariophyta_gene_name",
    "bacillariophyta_evalue", "bacillariophyta_bitscore", "bacillariophyta_pident",
    "bacillariophyta_qcov_percent", "bacillariophyta_scov_percent", "bacillariophyta_hit_confidence",
    "interproscan_analyses", "signature_accessions", "signature_descriptions", "interpro_accessions",
    "interpro_descriptions", "go_terms", "pathway_annotations", "transcriptome_descriptions",
    "transcriptome_broad_categories", "contig_id", "gene_model_start", "gene_model_end", "strand"
]
keep_cols = [c for c in keep_cols if c in df.columns]
clean = df[keep_cols].copy().rename(columns={"braker_protein_id": "gene_id"})
clean.to_csv(outfile, sep="\t", index=False)

print("Input rows:", len(df))
print("Output rows:", len(clean))
print("Output columns:", len(clean.columns))
print("Phaeodactylum yes/no counts:")
print(clean["in_Phaeodactylum_tricornutum"].value_counts(dropna=False))
print("Output file:", outfile)
