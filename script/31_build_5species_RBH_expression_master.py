#!/usr/bin/env python3
"""Build the Deer Lake expression + NI/SR/PT/TP RBH master table."""

import argparse
import csv

SPECIES = ("NI", "SR", "PT", "TP")

def read_table(path, key):
    out = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row[key]] = row
    return out

def choose(row, names):
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return ""

def read_map(path):
    out = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein = choose(row, [
                "protein", "protein_id", "Protein", "Protein_ID",
                "accession", "protein_accession"
            ])
            if protein:
                out[protein] = row
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--expression", required=True)
    for sp in SPECIES:
        p.add_argument(f"--{sp.lower()}-rbh", required=True)
        p.add_argument(f"--{sp.lower()}-map", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    expression = {}
    with open(args.expression) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            expression[row["DL_ORF"]] = row["DL_Average_TPM"]

    rbhs = {}
    maps = {}
    for sp in SPECIES:
        rbhs[sp] = read_table(getattr(args, f"{sp.lower()}_rbh"), "DL_ORF")
        maps[sp] = read_map(getattr(args, f"{sp.lower()}_map"))

    fields = ["DL_ORF", "DL_Average_TPM"]
    for sp in SPECIES:
        fields.extend([
            f"{sp}_protein", f"{sp}_transcript", f"{sp}_gene",
            f"{sp}_pident", f"{sp}_DL_query_coverage",
            f"{sp}_evalue", f"{sp}_bitscore"
        ])
    fields.extend([
        "NI_RBH", "SR_RBH", "PT_RBH", "TP_RBH",
        "RBH_pattern", "RBH_count"
    ])

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for orf, tpm in expression.items():
            row = {"DL_ORF": orf, "DL_Average_TPM": tpm}
            pattern = []

            for sp in SPECIES:
                hit = rbhs[sp].get(orf)
                present = 1 if hit else 0
                pattern.append(str(present))
                row[f"{sp}_RBH"] = present

                if hit:
                    protein = hit["reference_protein"]
                    meta = maps[sp].get(protein, {})
                    row[f"{sp}_protein"] = protein
                    row[f"{sp}_transcript"] = choose(meta, [
                        "transcript", "transcript_id", "Transcript", "Transcript_ID"
                    ])
                    row[f"{sp}_gene"] = choose(meta, [
                        "gene", "gene_id", "locus_tag", "Gene", "Gene_ID"
                    ])
                    row[f"{sp}_pident"] = hit.get("pident", "")
                    row[f"{sp}_DL_query_coverage"] = hit.get("DL_query_coverage", "")
                    row[f"{sp}_evalue"] = hit.get("evalue", "")
                    row[f"{sp}_bitscore"] = hit.get("bitscore", "")
                else:
                    for suffix in (
                        "protein", "transcript", "gene", "pident",
                        "DL_query_coverage", "evalue", "bitscore"
                    ):
                        row[f"{sp}_{suffix}"] = ""

            row["RBH_pattern"] = "".join(pattern)
            row["RBH_count"] = sum(int(x) for x in pattern)
            writer.writerow(row)

    print("Deer Lake ORFs:", len(expression))
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
