#!/usr/bin/env python3
"""Add all Phaeodactylum tricornutum BLASTN gene-overlap hits to the BRAKER4 table."""

import pandas as pd
import re

master_table = "DL_diatom_BRAKER_annotation_expression_lengths.tsv"
diatom_pt_hits_file = "/work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast_18_diatom_v1/blast_out/18_diatom_nuclear_v1_genes_with_phaeodactylum_hits.tsv"
pt_gene_map_file = "/work/ebg_lab/eb/diatom_consortia/phaeodactylum_blast_18_diatom_v1/blast_out/blast_hit_to_phaeodactylum_gene.tsv"
out_file = "DL_diatom_BRAKER_annotation_expression_lengths_ALL_Phaeodactylum_hits.tsv"

def gene_root(x):
    if pd.isna(x):
        return ""
    return re.sub(r"\.t[0-9]+$", "", str(x).strip())

def join_unique(x):
    return ";".join(sorted(set(x.dropna().astype(str))))

master = pd.read_csv(master_table, sep="\t", low_memory=False)
master.columns = master.columns.str.strip()
master["braker_gene_root"] = master["braker_protein_id"].apply(gene_root)

hit_cols = [
    "diatom_contig", "diatom_gene_start", "diatom_gene_end", "diatom_gene_id",
    "diatom_gene_score", "diatom_gene_strand", "hit_contig", "hit_start", "hit_end",
    "blast_hit_id", "bitscore", "blast_strand", "pt_contig", "pt_start", "pt_end",
    "pident", "aln_len", "evalue"
]
hits = pd.read_csv(diatom_pt_hits_file, sep="\t", names=hit_cols, low_memory=False)
hits["diatom_gene_root_for_merge"] = hits["diatom_gene_id"].apply(gene_root)

pt_map = pd.read_csv(pt_gene_map_file, sep="\t", low_memory=False)
pt_map.columns = pt_map.columns.str.strip()
hits = hits.merge(pt_map, on="blast_hit_id", how="left")

pt_summary = (
    hits.groupby("diatom_gene_root_for_merge", as_index=False)
    .agg(
        phaeodactylum_num_all_hits=("blast_hit_id", "count"),
        phaeodactylum_num_unique_hit_ids=("blast_hit_id", "nunique"),
        phaeodactylum_num_unique_PT_genes=("pt_gene_id", "nunique"),
        phaeodactylum_all_blast_hit_ids=("blast_hit_id", join_unique),
        phaeodactylum_all_pt_gene_ids=("pt_gene_id", join_unique),
        phaeodactylum_all_pt_gene_names=("pt_gene_name", join_unique),
        phaeodactylum_all_pt_locus_tags=("pt_locus_tag", join_unique),
        phaeodactylum_best_bitscore=("bitscore", "max"),
        phaeodactylum_best_pident=("pident", "max"),
        phaeodactylum_best_evalue=("evalue", "min"),
        phaeodactylum_max_alignment_length=("aln_len", "max")
    )
)

merged = master.merge(pt_summary, left_on="braker_gene_root", right_on="diatom_gene_root_for_merge", how="left")
merged["in_Phaeodactylum_tricornutum"] = merged["phaeodactylum_num_all_hits"].notna().map({True: "yes", False: "no"})
merged = merged.drop(columns=["diatom_gene_root_for_merge"], errors="ignore")

front_cols = [
    "braker_protein_id", "compartment", "gene_model_length_bp", "cds_length_bp", "protein_length_aa",
    "recommended_annotation", "recommended_annotation_source", "recommended_annotation_confidence",
    "average_expression_TPM", "in_Phaeodactylum_tricornutum", "phaeodactylum_num_all_hits",
    "phaeodactylum_num_unique_PT_genes", "phaeodactylum_all_pt_gene_ids", "phaeodactylum_all_pt_gene_names",
    "phaeodactylum_all_pt_locus_tags", "phaeodactylum_best_bitscore", "phaeodactylum_best_pident",
    "phaeodactylum_best_evalue", "phaeodactylum_max_alignment_length", "has_transcriptome_ORF_hit",
    "num_unique_ORFs", "mapped_orf_ids", "has_any_annotation", "has_swissprot_hit",
    "has_bacillariophyta_hit", "has_interproscan_hit"
]
front_cols = [c for c in front_cols if c in merged.columns]
merged = merged[front_cols + [c for c in merged.columns if c not in front_cols]]
merged.to_csv(out_file, sep="\t", index=False)

print("Input master rows:", len(master))
print("All Phaeodactylum hit rows:", len(hits))
print("Unique Deer Lake genes with Phaeodactylum hits:", pt_summary["diatom_gene_root_for_merge"].nunique())
print("Output rows:", len(merged))
print("Genes/protein isoforms marked yes:", merged["in_Phaeodactylum_tricornutum"].eq("yes").sum())
print("Genes/protein isoforms marked no:", merged["in_Phaeodactylum_tricornutum"].eq("no").sum())
print("Output file:", out_file)
