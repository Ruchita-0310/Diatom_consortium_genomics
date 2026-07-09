#!/usr/bin/env python3
"""
Create the final simple read-level table for mixed diatom-bacterial Hi-C pairs.

Each mixed Hi-C read pair is represented by two rows:
  - read_number = 1 for HiC_R1
  - read_number = 2 for HiC_R2

Run from:
  /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
"""

from pathlib import Path
import pandas as pd

INFILE = Path("03_tables/HiC_DIATOM_BACTERIAL_read_level_MAPQ30_PID95.tsv")
OUTFILE = Path("03_tables/HiC_DIATOM_BACTERIAL_read_level_SIMPLE_COLUMNS_MAPQ30_PID95.tsv")


def main() -> None:
    df = pd.read_csv(INFILE, sep="\t")

    read1 = pd.DataFrame({
        "read_id": df["read_id"],
        "read_number": 1,
        "contig_id": df["forward_contig_id"],
        "paired_contig_id": df["reverse_contig_id"],
        "mapq": df["forward_mapq"],
        "percent_identity": df["forward_percent_identity"],
        "aligned_length_bp": df["forward_aligned_length_bp"],
        "pair_type_code": df["pair_type_code"],
        "pair_type": df["pair_type"],
    })

    read2 = pd.DataFrame({
        "read_id": df["read_id"],
        "read_number": 2,
        "contig_id": df["reverse_contig_id"],
        "paired_contig_id": df["forward_contig_id"],
        "mapq": df["reverse_mapq"],
        "percent_identity": df["reverse_percent_identity"],
        "aligned_length_bp": df["reverse_aligned_length_bp"],
        "pair_type_code": df["pair_type_code"],
        "pair_type": df["pair_type"],
    })

    simple = pd.concat([read1, read2], ignore_index=True).sort_values(["read_id", "read_number"])
    simple.to_csv(OUTFILE, sep="\t", index=False)

    unique_pairs = simple["read_id"].nunique()
    print("Created:", OUTFILE)
    print("Unique mixed Hi-C read pairs:", unique_pairs)
    print("Read rows:", len(simple))
    print("Lines including header:", len(simple) + 1)
    print(simple.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
