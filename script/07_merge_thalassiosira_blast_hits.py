#!/usr/bin/env python3
"""Merge raw Thalassiosira-vs-diatom BLASTN hits with gene overlaps.

This mirrors the Phaeodactylum comparison workflow:
  1. one row per raw BLASTN hit;
  2. overlapping reference genes collapsed with semicolons;
  3. overlapping Deer Lake BRAKER4 genes collapsed with semicolons;
  4. no BLASTN filtering beyond the parameters used when BLASTN was run.
"""

from __future__ import annotations

from pathlib import Path
import sys
import pandas as pd

RAW_BLAST = Path("02_blast/thalassiosira_vs_diatom_dcmegablast.tsv")
THAL_OVERLAP = Path("03_filtered/thalassiosira_hits_with_Thalassiosira_genes.tsv")
DIATOM_OVERLAP = Path("03_filtered/diatom_hits_with_BRAKER_ET_genes.tsv")
OUTPUT = Path("04_summary/thalassiosira_vs_diatom_BLASTN_with_Thalassiosira_and_BRAKER_ET_genes.tsv")
SUMMARY = Path("04_summary/thalassiosira_vs_diatom_BLASTN_summary.txt")

BLAST_COLUMNS = [
    "thal_contig", "diatom_contig", "pident", "aln_len", "mismatch",
    "gapopen", "thal_start", "thal_end", "diatom_start", "diatom_end",
    "evalue", "bitscore", "thal_len", "diatom_len", "qcovs",
]

THAL_INTERSECT_COLUMNS = [
    "a_thal_contig", "a_start", "a_end", "blast_hit_id", "a_pident",
    "a_aln_len", "a_evalue", "a_bitscore", "a_thal_len",
    "a_diatom_contig", "a_diatom_start", "a_diatom_end", "a_hit_strand",
    "b_thal_contig", "b_gene_start", "b_gene_end", "thal_gene_id",
    "thal_gene_name", "thal_gene_symbol", "thal_locus_tag",
    "thal_gene_strand",
]

DIATOM_INTERSECT_COLUMNS = [
    "a_diatom_contig", "a_start", "a_end", "blast_hit_id", "a_pident",
    "a_aln_len", "a_evalue", "a_bitscore", "a_diatom_len",
    "a_thal_contig", "a_thal_start", "a_thal_end", "a_hit_strand",
    "b_diatom_contig", "b_gene_start", "b_gene_end", "diatom_gene_id",
    "diatom_gene_attr_id", "diatom_gene_strand",
]

MISSING_TOKENS = {"", ".", "-1", "NA", "NaN", "nan", "None"}


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required input file not found: {path}")


def clean_value(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text in MISSING_TOKENS:
        return None
    return text


def collapse_unique(values: pd.Series) -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for raw in values:
        value = clean_value(raw)
        if value is None or value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ";".join(ordered) if ordered else "NA"


def count_unique_semicolon(series: pd.Series) -> int:
    values: set[str] = set()
    for raw in series:
        value = clean_value(raw)
        if value is None:
            continue
        for item in value.split(";"):
            item = item.strip()
            if item and item not in MISSING_TOKENS:
                values.add(item)
    return len(values)


def main() -> int:
    for path in (RAW_BLAST, THAL_OVERLAP, DIATOM_OVERLAP):
        require_file(path)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    blast = pd.read_csv(RAW_BLAST, sep="\t", header=None, names=BLAST_COLUMNS, dtype=str)
    blast.insert(0, "blast_hit_id", [f"hit_{i:06d}" for i in range(1, len(blast) + 1)])

    thal = pd.read_csv(
        THAL_OVERLAP, sep="\t", header=None, names=THAL_INTERSECT_COLUMNS,
        dtype=str, keep_default_na=False,
    )
    diatom = pd.read_csv(
        DIATOM_OVERLAP, sep="\t", header=None, names=DIATOM_INTERSECT_COLUMNS,
        dtype=str, keep_default_na=False,
    )

    thal_fields = [
        "thal_gene_id", "thal_gene_name", "thal_gene_symbol",
        "thal_locus_tag", "thal_gene_strand",
    ]
    diatom_fields = [
        "diatom_gene_id", "diatom_gene_attr_id", "diatom_gene_strand",
    ]

    thal_grouped = (
        thal.groupby("blast_hit_id", sort=False)[thal_fields]
        .agg(collapse_unique)
        .reset_index()
    )
    diatom_grouped = (
        diatom.groupby("blast_hit_id", sort=False)[diatom_fields]
        .agg(collapse_unique)
        .reset_index()
    )

    merged = blast.merge(thal_grouped, on="blast_hit_id", how="left")
    merged = merged.merge(diatom_grouped, on="blast_hit_id", how="left")

    annotation_fields = thal_fields + diatom_fields
    merged[annotation_fields] = merged[annotation_fields].fillna("NA")

    if len(merged) != len(blast):
        raise RuntimeError(
            f"Expected one output row per BLAST hit ({len(blast)}), got {len(merged)}"
        )

    merged.to_csv(OUTPUT, sep="\t", index=False)

    thal_present = merged["thal_gene_id"].ne("NA")
    diatom_present = merged["diatom_gene_id"].ne("NA")
    both_present = thal_present & diatom_present

    summary_lines = [
        f"total_blast_hits\t{len(merged):,}",
        f"hits_with_Thalassiosira_gene\t{int(thal_present.sum()):,}",
        f"hits_with_diatom_BRAKER_ET_gene\t{int(diatom_present.sum()):,}",
        f"hits_with_both_Thalassiosira_and_diatom_gene\t{int(both_present.sum()):,}",
        f"unique_Thalassiosira_genes_hit\t{count_unique_semicolon(merged.loc[thal_present, 'thal_gene_id']):,}",
        f"unique_diatom_BRAKER_ET_genes_hit\t{count_unique_semicolon(merged.loc[diatom_present, 'diatom_gene_id']):,}",
    ]
    SUMMARY.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"Wrote: {OUTPUT}")
    print(f"Rows: {len(merged):,}")
    print(f"Columns: {len(merged.columns)}")
    print(f"Wrote: {SUMMARY}")
    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
