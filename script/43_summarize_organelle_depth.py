#!/usr/bin/env python3
"""Summarize per-base samtools depth for organelle candidate contigs."""

import argparse
import statistics
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize a samtools depth -aa table by contig. "
            "Reports length, mean depth, median depth, minimum and maximum "
            "depth, covered bases, and breadth of coverage."
        )
    )
    parser.add_argument("depth_tsv", help="samtools depth -aa output")
    return parser.parse_args()


def main():
    args = parse_args()
    depths = defaultdict(list)

    with open(args.depth_tsv) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Unexpected depth line: {line.rstrip()}")
            contig = fields[0]
            depth = int(fields[2])
            depths[contig].append(depth)

    print(
        "Contig\tLength_bp\tMean_depth\tMedian_depth\tMin_depth\t"
        "Max_depth\tCovered_bases\tBreadth_percent"
    )

    for contig, values in depths.items():
        length = len(values)
        covered = sum(d > 0 for d in values)
        breadth = (covered / length * 100.0) if length else 0.0
        mean_depth = statistics.fmean(values) if values else 0.0
        median_depth = statistics.median(values) if values else 0.0

        print(
            f"{contig}\t{length}\t{mean_depth:.2f}\t{median_depth:g}\t"
            f"{min(values)}\t{max(values)}\t{covered}\t{breadth:.2f}"
        )


if __name__ == "__main__":
    main()
