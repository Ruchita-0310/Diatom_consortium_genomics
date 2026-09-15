#!/usr/bin/env python3
"""Prepared next step: create a 30-gene table with five representatives per system.

Verify organelle versus nuclear origin of NODE_6033.p2 before final figure use.
"""

import argparse
import csv

SELECTED = [
    "NODE_67.p4", "NODE_63.p1", "NODE_44.p4", "NODE_30207.p1", "NODE_7928.p2",
    "NODE_6033.p2", "NODE_8468.p1", "NODE_18074.p1", "NODE_20723.p1", "NODE_10495.p2",
    "NODE_8343.p1", "NODE_5705.p1", "NODE_9414.p1", "NODE_10595.p1", "NODE_10596.p1",
    "NODE_14611.p3", "NODE_5391.p2", "NODE_5837.p1", "NODE_1090.p2", "NODE_4370.p1",
    "NODE_971.p2", "NODE_8419.p1", "NODE_12386.p1", "NODE_21007.p2", "NODE_10703.p1",
    "NODE_2665.p1", "NODE_10018.p1", "NODE_25198.p1", "NODE_23676.p1", "NODE_26832.p1",
]

ORGANELLE_NOT_ASSESSED = {"NODE_67.p4", "NODE_63.p1", "NODE_44.p4"}

FIELDS = [
    "Functional_category", "Curated_function", "DL_ORF",
    "DL_Average_TPM", "DL_expression_percentile",
    "NI_RBH", "SR_RBH", "PT_RBH", "TP_RBH",
    "RBH_pattern", "RBH_count", "RBH_assessment",
    "KOfam_KO", "KOfam_definition",
    "EggNOG_preferred_name", "EggNOG_description", "Direct_Pfam",
]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = {}
    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows[row["DL_ORF"]] = row

    missing = [orf for orf in SELECTED if orf not in rows]
    if missing:
        raise SystemExit("Missing selected ORFs: " + ", ".join(missing))

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()

        for orf in SELECTED:
            row = rows[orf].copy()

            if orf in ORGANELLE_NOT_ASSESSED:
                for column in (
                    "NI_RBH", "SR_RBH", "PT_RBH", "TP_RBH",
                    "RBH_pattern", "RBH_count"
                ):
                    row[column] = "NA"
                row["RBH_assessment"] = "Not assessed: plastid encoded gene"
            else:
                row["RBH_assessment"] = "Reciprocal best protein hit"

            writer.writerow(row)

    print("Selected genes:", len(SELECTED))
    print("Expected:", 30)
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
