#!/usr/bin/env python3
"""
Parse separate Hi-C R1/R2 BAM files and keep high-confidence primary alignments.

Filters:
  - primary alignments only: removes unmapped, secondary, and supplementary alignments
  - MAPQ >= 30
  - percent identity >= 95, calculated from NM tag and aligned CIGAR length

Run from:
  /work/ebg_lab/eb/diatom_consortia/hic_bwa_separate_reads
"""

import re
import subprocess
from pathlib import Path

MIN_MAPQ = 30
MIN_PIDENT = 95.0
OUTDIR = Path("03_tables")
OUTDIR.mkdir(exist_ok=True)


def aligned_length_from_cigar(cigar: str) -> int:
    """Return aligned length from CIGAR operations M, =, and X."""
    total = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in ("M", "=", "X"):
            total += int(n)
    return total


def make_table(bam: str, out_tsv: str) -> None:
    kept = 0
    total = 0

    cmd = [
        "samtools", "view",
        "-F", "2308",      # remove unmapped, secondary, supplementary
        "-q", str(MIN_MAPQ),
        bam,
    ]

    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc, open(out_tsv, "w") as out:
        out.write("read_id\tcontig_id\tstart\tmapq\tcigar\tNM\taligned_length_bp\tpercent_identity\n")

        for line in proc.stdout:
            total += 1
            fields = line.rstrip("\n").split("\t")
            read_id = fields[0]
            contig_id = fields[2]
            start = fields[3]
            mapq = fields[4]
            cigar = fields[5]

            nm = None
            for tag in fields[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag.split(":")[-1])
                    break
            if nm is None:
                continue

            aln_len = aligned_length_from_cigar(cigar)
            if aln_len == 0:
                continue

            pident = 100.0 * (aln_len - nm) / aln_len
            if pident >= MIN_PIDENT:
                out.write(
                    f"{read_id}\t{contig_id}\t{start}\t{mapq}\t{cigar}\t{nm}\t{aln_len}\t{pident:.3f}\n"
                )
                kept += 1

    print(f"{out_tsv}: kept {kept:,} high-identity reads from {total:,} primary MAPQ>={MIN_MAPQ} alignments")


def main() -> None:
    make_table("02_alignments/HiC_R1.sorted.bam", "03_tables/HiC_R1.primary_MAPQ30_PID95.tsv")
    make_table("02_alignments/HiC_R2.sorted.bam", "03_tables/HiC_R2.primary_MAPQ30_PID95.tsv")


if __name__ == "__main__":
    main()
