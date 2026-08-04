#!/usr/bin/env python3
"""
Join high-confidence Hi-C R1 and R2 alignments by read ID and classify pair type.

Contig classification rule:
  - contigs present in 18_diatom.fasta = diatom
  - all other contigs in the polished whole assembly = bacterial

Run from:
  /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
"""

from pathlib import Path
import pandas as pd

DIATOM_FASTA = Path("/work/ebg_lab/eb/diatom_consortia/metatranscriptomics/genome_index/18_diatom.fasta")
WHOLE_FASTA = Path("/work/ebg_lab/eb/diatom_consortia/MAGS_guppy/1_sr_pypolca_output/pypolca_corrected.fasta")

R1_FILE = Path("03_tables/HiC_R1.primary_MAPQ30_PID95.tsv")
R2_FILE = Path("03_tables/HiC_R2.primary_MAPQ30_PID95.tsv")
OUTDIR = Path("03_tables")


def fasta_ids(path: Path) -> set[str]:
    ids = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def make_contig_type_map() -> pd.DataFrame:
    diatom_ids = fasta_ids(DIATOM_FASTA)
    rows = []
    with open(WHOLE_FASTA) as handle:
        for line in handle:
            if line.startswith(">"):
                contig = line[1:].split()[0]
                rows.append({
                    "contig_id": contig,
                    "contig_type": "diatom" if contig in diatom_ids else "bacterial",
                })
    contig_map = pd.DataFrame(rows)
    contig_map.to_csv(OUTDIR / "whole_assembly_contig_type_map.tsv", sep="\t", index=False)
    return contig_map


def classify_pair(row: pd.Series) -> tuple[int, str]:
    f = row["forward_contig_type"]
    r = row["reverse_contig_type"]

    if f == "diatom" and r == "diatom":
        return 1, "both_diatom"
    if f == "bacterial" and r == "bacterial":
        return 2, "both_bacterial"
    if f == "diatom" and r == "bacterial":
        return 3, "forward_diatom_reverse_bacterial"
    if f == "bacterial" and r == "diatom":
        return 4, "forward_bacterial_reverse_diatom"
    return 9, "unknown"


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)

    contig_map = make_contig_type_map()
    type_dict = dict(zip(contig_map["contig_id"], contig_map["contig_type"]))

    r1 = pd.read_csv(R1_FILE, sep="\t").rename(columns={
        "contig_id": "forward_contig_id",
        "start": "forward_start",
        "mapq": "forward_mapq",
        "cigar": "forward_cigar",
        "NM": "forward_NM",
        "aligned_length_bp": "forward_aligned_length_bp",
        "percent_identity": "forward_percent_identity",
    })

    r2 = pd.read_csv(R2_FILE, sep="\t").rename(columns={
        "contig_id": "reverse_contig_id",
        "start": "reverse_start",
        "mapq": "reverse_mapq",
        "cigar": "reverse_cigar",
        "NM": "reverse_NM",
        "aligned_length_bp": "reverse_aligned_length_bp",
        "percent_identity": "reverse_percent_identity",
    })

    pairs = r1.merge(r2, on="read_id", how="inner")
    pairs["forward_contig_type"] = pairs["forward_contig_id"].map(type_dict).fillna("unknown")
    pairs["reverse_contig_type"] = pairs["reverse_contig_id"].map(type_dict).fillna("unknown")

    classified = pairs.apply(classify_pair, axis=1, result_type="expand")
    pairs["pair_type_code"] = classified[0]
    pairs["pair_type"] = classified[1]

    pairs.to_csv(OUTDIR / "HiC_read_pairs_MAPQ30_PID95_with_contig_types.tsv", sep="\t", index=False)

    summary = (
        pairs.groupby(["pair_type_code", "pair_type"], as_index=False)
        .size()
        .rename(columns={"size": "read_pair_count"})
        .sort_values("pair_type_code")
    )
    summary["percent"] = 100 * summary["read_pair_count"] / summary["read_pair_count"].sum()
    summary.to_csv(OUTDIR / "HiC_pair_type_summary_MAPQ30_PID95.tsv", sep="\t", index=False)

    mixed = pairs[pairs["pair_type_code"].isin([3, 4])].copy()
    mixed["diatom_contig_id"] = mixed.apply(
        lambda row: row["forward_contig_id"] if row["forward_contig_type"] == "diatom" else row["reverse_contig_id"],
        axis=1,
    )
    mixed["bacterial_contig_id"] = mixed.apply(
        lambda row: row["forward_contig_id"] if row["forward_contig_type"] == "bacterial" else row["reverse_contig_id"],
        axis=1,
    )
    mixed.to_csv(OUTDIR / "HiC_DIATOM_BACTERIAL_read_level_MAPQ30_PID95.tsv", sep="\t", index=False)

    print("R1 high-identity reads:", len(r1))
    print("R2 high-identity reads:", len(r2))
    print("Joined read pairs:", len(pairs))
    print("Mixed diatom-bacterial read pairs:", len(mixed))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
