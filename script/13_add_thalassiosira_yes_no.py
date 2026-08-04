#!/usr/bin/env python3
"""Add a Thalassiosira pseudonana gene-linked BLASTN yes/no field.

The script preserves all existing columns, including the prior
present_in_Phaeodactylum_tricornutum field, and adds:
    present_in_Thalassiosira_pseudonana

A Deer Lake BRAKER4 gene is marked yes only when the same raw BLASTN hit
is linked to both an annotated T. pseudonana gene and a Deer Lake BRAKER4 gene.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
import pandas as pd

PRESENT_COLUMN = "present_in_Thalassiosira_pseudonana"
MISSING = {"", ".", "NA", "NaN", "nan", "None"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-final",
        default="09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table.tsv",
        help="Existing clean BRAKER4 isoform-level table.",
    )
    parser.add_argument(
        "--blast-summary",
        default=(
            "/work/ebg_lab/eb/diatom_consortia/"
            "thalassiosira_to_diatom_blastn_redo/04_summary/"
            "thalassiosira_vs_diatom_BLASTN_with_Thalassiosira_and_BRAKER_ET_genes.tsv"
        ),
        help="Gene-linked Thalassiosira-vs-diatom BLASTN table.",
    )
    parser.add_argument(
        "--output-final",
        default="09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP.tsv",
    )
    parser.add_argument(
        "--output-sorted",
        default="09_final/DL_diatom_FINAL_clean_BRAKER_isoform_table_PT_TP_sorted_by_Average_TPM.tsv",
    )
    parser.add_argument(
        "--output-review",
        default="09_final/DL_diatom_FINAL_gene_table_for_boss_PT_TP.tsv",
    )
    return parser.parse_args()


def find_column(df: pd.DataFrame, candidates: list[str], description: str) -> str:
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    raise KeyError(
        f"Could not identify {description}. Tried: {', '.join(candidates)}. "
        f"Available columns: {', '.join(df.columns)}"
    )


def clean_gene_id(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text in MISSING:
        return None
    text = re.sub(r"^gene-", "", text)
    return text


def gene_root(value: object) -> str | None:
    text = clean_gene_id(value)
    if text is None:
        return None
    # BRAKER isoforms are typically g1234.t1; the comparison table stores g1234.
    return re.sub(r"\.t\d+$", "", text)


def split_gene_field(value: object) -> list[str]:
    if pd.isna(value):
        return []
    output: list[str] = []
    for item in str(value).split(";"):
        cleaned = clean_gene_id(item)
        if cleaned is not None:
            output.append(cleaned)
    return output


def main() -> int:
    args = parse_args()
    input_final = Path(args.input_final)
    blast_summary = Path(args.blast_summary)
    output_final = Path(args.output_final)
    output_sorted = Path(args.output_sorted)
    output_review = Path(args.output_review)

    for path in (input_final, blast_summary):
        if not path.is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    final = pd.read_csv(input_final, sep="\t", dtype=str, keep_default_na=False)
    blast = pd.read_csv(blast_summary, sep="\t", dtype=str, keep_default_na=False)

    for required in ("thal_gene_id", "diatom_gene_id"):
        if required not in blast.columns:
            raise KeyError(f"Missing required BLAST summary column: {required}")

    linked_roots: set[str] = set()
    for _, row in blast.iterrows():
        thal_genes = split_gene_field(row["thal_gene_id"])
        diatom_genes = split_gene_field(row["diatom_gene_id"])
        if not thal_genes or not diatom_genes:
            continue
        for gene in diatom_genes:
            root = gene_root(gene)
            if root is not None:
                linked_roots.add(root)

    id_column = find_column(final, ["protein_id", "gene_id", "transcript_id"], "BRAKER isoform ID column")
    final[PRESENT_COLUMN] = final[id_column].map(
        lambda value: "yes" if gene_root(value) in linked_roots else "no"
    )

    output_final.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(output_final, sep="\t", index=False)

    tpm_column = next(
        (column for column in ["Average_TPM", "diatom_Average_TPM"] if column in final.columns),
        None,
    )
    if tpm_column is not None:
        sorted_final = final.copy()
        sorted_final["__TPM_numeric"] = pd.to_numeric(sorted_final[tpm_column], errors="coerce")
        sorted_final = sorted_final.sort_values(
            by="__TPM_numeric", ascending=False, na_position="last"
        ).drop(columns="__TPM_numeric")
        sorted_final.to_csv(output_sorted, sep="\t", index=False)
    else:
        final.to_csv(output_sorted, sep="\t", index=False)
        print("WARNING: Average_TPM column not found; sorted output preserves input order.")

    review_mapping = {
        "gene_id": id_column,
        "contig_id": find_column(final, ["contig_id", "seqid"], "contig column"),
        "diatom_compartment": find_column(final, ["diatom_compartment", "compartment"], "compartment column"),
        "diatom_gene_length_bp": find_column(final, ["diatom_gene_length_bp", "gene_length_bp", "gene_length"], "gene length column"),
        "functional_annotation": find_column(final, ["functional_annotation", "recommended_annotation"], "functional annotation column"),
        "diatom_Average_TPM": find_column(final, ["diatom_Average_TPM", "Average_TPM"], "Average_TPM column"),
    }

    review = pd.DataFrame({new: final[old] for new, old in review_mapping.items()})
    pt_column = "present_in_Phaeodactylum_tricornutum"
    if pt_column in final.columns:
        review[pt_column] = final[pt_column]
    review[PRESENT_COLUMN] = final[PRESENT_COLUMN]
    review.to_csv(output_review, sep="\t", index=False)

    counts = final[PRESENT_COLUMN].value_counts(dropna=False)
    print(f"Input rows: {len(final):,}")
    print(f"Unique Deer Lake gene roots linked to annotated T. pseudonana genes: {len(linked_roots):,}")
    print(f"Wrote: {output_final}")
    print(f"Wrote: {output_sorted}")
    print(f"Wrote: {output_review}")
    print("Thalassiosira yes/no counts:")
    for label, count in counts.items():
        print(f"  {label}: {count:,}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
