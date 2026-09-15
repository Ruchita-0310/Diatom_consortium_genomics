#!/usr/bin/env python3
"""Extract reciprocal best BLASTP hits using unique maximum bitscore."""

import argparse
import csv
from collections import defaultdict

FIELDS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

def unique_top_hits(path):
    by_query = defaultdict(list)
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=FIELDS)
        for row in reader:
            row["bitscore_float"] = float(row["bitscore"])
            by_query[row["qseqid"]].append(row)

    best = {}
    for query, rows in by_query.items():
        top_score = max(r["bitscore_float"] for r in rows)
        top_rows = [r for r in rows if r["bitscore_float"] == top_score]
        if len({r["sseqid"] for r in top_rows}) == 1:
            best[query] = top_rows[0]
    return best

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--forward", required=True)
    p.add_argument("--reverse", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    forward = unique_top_hits(args.forward)
    reverse = unique_top_hits(args.reverse)

    fields = [
        "DL_ORF", "reference_protein", "pident",
        "DL_query_coverage", "evalue", "bitscore"
    ]

    n = 0
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for dl_orf, hit in forward.items():
            ref = hit["sseqid"]
            rev = reverse.get(ref)
            if rev is None or rev["sseqid"] != dl_orf:
                continue

            qlen = float(hit["qlen"])
            aln_len = float(hit["length"])
            qcov = 100.0 * aln_len / qlen if qlen else 0.0

            writer.writerow({
                "DL_ORF": dl_orf,
                "reference_protein": ref,
                "pident": hit["pident"],
                "DL_query_coverage": f"{qcov:.3f}",
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
            })
            n += 1

    print("Unique reciprocal best-hit pairs:", n)
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
