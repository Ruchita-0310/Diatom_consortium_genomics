#!/usr/bin/env python3
"""Integrate KOfam, EggNOG and direct Pfam/HMMER with DL expression + RBH."""

import argparse
import csv
import gzip
from collections import defaultdict

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--master", required=True)
    p.add_argument("--eggnog", required=True)
    p.add_argument("--kofam", required=True)
    p.add_argument("--pfam", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    eggnog = {}
    with gzip.open(args.eggnog, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            eggnog[row["orf"]] = {
                "EggNOG_description": row.get("description", ""),
                "EggNOG_preferred_name": row.get("preferred_name", ""),
                "EggNOG_GO": row.get("gos", ""),
                "EggNOG_EC": row.get("ec", ""),
                "EggNOG_KEGG_KO": row.get("kegg_ko", ""),
                "EggNOG_KEGG_pathway": row.get("kegg_pathway", ""),
                "EggNOG_Pfam": row.get("pfams", ""),
            }

    kofam = {}
    with gzip.open(args.kofam, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            kofam[row["orf"]] = {
                "KOfam_KO": row.get("ko", ""),
                "KOfam_definition": row.get("ko_definition", ""),
                "KOfam_threshold": row.get("thrshld", ""),
                "KOfam_score": row.get("score", ""),
                "KOfam_evalue": row.get("evalue", ""),
            }

    pfam_hits = defaultdict(list)
    with gzip.open(args.pfam, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split(maxsplit=18)
            if len(fields) < 18:
                continue

            orf = fields[0]
            pfam_name = fields[2]
            pfam_accession = fields[3]

            try:
                included_domains = int(fields[17])
            except ValueError:
                continue

            if included_domains <= 0:
                continue

            hit = f"{pfam_accession}|{pfam_name}"
            if hit not in pfam_hits[orf]:
                pfam_hits[orf].append(hit)

    with open(args.master) as inp, open(args.output, "w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")

        new_columns = [
            "KOfam_KO", "KOfam_definition", "KOfam_threshold",
            "KOfam_score", "KOfam_evalue",
            "EggNOG_preferred_name", "EggNOG_description",
            "EggNOG_KEGG_KO", "EggNOG_KEGG_pathway",
            "EggNOG_EC", "EggNOG_GO", "EggNOG_Pfam",
            "Direct_Pfam",
        ]

        base = reader.fieldnames
        insert_after = base.index("DL_expression_percentile") + 1
        fields = base[:insert_after] + new_columns + base[insert_after:]

        writer = csv.DictWriter(
            out, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()

        total = 0
        any_annotation = 0

        for row in reader:
            total += 1
            orf = row["DL_ORF"]
            k = kofam.get(orf, {})
            e = eggnog.get(orf, {})

            row.update({
                "KOfam_KO": k.get("KOfam_KO", ""),
                "KOfam_definition": k.get("KOfam_definition", ""),
                "KOfam_threshold": k.get("KOfam_threshold", ""),
                "KOfam_score": k.get("KOfam_score", ""),
                "KOfam_evalue": k.get("KOfam_evalue", ""),
                "EggNOG_preferred_name": e.get("EggNOG_preferred_name", ""),
                "EggNOG_description": e.get("EggNOG_description", ""),
                "EggNOG_KEGG_KO": e.get("EggNOG_KEGG_KO", ""),
                "EggNOG_KEGG_pathway": e.get("EggNOG_KEGG_pathway", ""),
                "EggNOG_EC": e.get("EggNOG_EC", ""),
                "EggNOG_GO": e.get("EggNOG_GO", ""),
                "EggNOG_Pfam": e.get("EggNOG_Pfam", ""),
                "Direct_Pfam": ";".join(pfam_hits.get(orf, [])),
            })

            if row["KOfam_KO"] or row["EggNOG_description"] or row["Direct_Pfam"]:
                any_annotation += 1

            writer.writerow(row)

    print("EggNOG annotated ORFs:", len(eggnog))
    print("KOfam annotated ORFs:", len(kofam))
    print("ORFs with included Pfam hits:", len(pfam_hits))
    print("Total DL ORFs:", total)
    print("ORFs with at least one annotation source:", any_annotation)
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
