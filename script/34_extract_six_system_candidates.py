#!/usr/bin/env python3
"""Broad first-pass collection of candidates in six functional systems."""

import argparse
import csv
import re
from collections import defaultdict

PATTERNS = {
    "Photosynthesis": [
        r"photosystem", r"psa[a-z0-9]*\b", r"psb[a-z0-9]*\b",
        r"light[- ]harvest", r"chlorophyll", r"fucoxanthin",
        r"cytochrome b6", r"pet[a-z]\b", r"ferredoxin",
        r"ferredoxin.nadp", r"plastocyanin",
    ],
    "Carbon fixation & CCM": [
        r"rubisco", r"ribulose.*bisphosphate", r"carbonic anhydrase",
        r"bicarbonate", r"carbonate", r"carbon concentrating",
        r"\bccm\b", r"phosphoribulokinase",
        r"sedoheptulose.*bisphosphatase",
        r"fructose.*bisphosphatase",
        r"glyceraldehyde.*phosphate dehydrogenase",
    ],
    "Silica metabolism": [
        r"silicon transporter", r"silicate transporter", r"silica",
        r"silaffin", r"silacidin", r"frustulin", r"\bsit[123]?\b",
    ],
    "Ion homeostasis / osmoregulation": [
        r"sodium", r"potassium", r"proton", r"cation", r"anion",
        r"antiporter", r"exchanger", r"Na\+", r"K\+", r"Ca2\+",
        r"calcium", r"chloride", r"osmotic", r"osmoreg",
        r"vacuolar.*ATPase", r"V-type.*ATPase",
    ],
    "Urea & nitrogen metabolism": [
        r"urea", r"urease", r"urea transporter", r"ornithine",
        r"argininosuccinate", r"carbamoyl.*phosphate",
        r"nitrite", r"nitrate", r"ammonium", r"ammonia",
        r"glutamine synthetase", r"glutamate synthase",
        r"glutamate dehydrogenase", r"nitrogen metabolism",
    ],
    "Oxidative stress": [
        r"superoxide dismutase", r"catalase", r"peroxidase",
        r"peroxiredoxin", r"thioredoxin", r"glutaredoxin",
        r"glutathione", r"ascorbate", r"oxidative stress",
        r"reactive oxygen", r"\bros\b",
    ],
}

def annotation_text(row):
    return " | ".join([
        row.get("KOfam_KO", ""),
        row.get("KOfam_definition", ""),
        row.get("EggNOG_preferred_name", ""),
        row.get("EggNOG_description", ""),
        row.get("EggNOG_KEGG_KO", ""),
        row.get("EggNOG_KEGG_pathway", ""),
        row.get("EggNOG_Pfam", ""),
        row.get("Direct_Pfam", ""),
    ])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    hits = []
    counts = defaultdict(int)

    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            text = annotation_text(row)
            categories = []
            matched_terms = []

            for category, patterns in PATTERNS.items():
                local = [pat for pat in patterns if re.search(pat, text, re.I)]
                if local:
                    categories.append(category)
                    matched_terms.extend(local)

            if not categories:
                continue

            row["Functional_category_candidate"] = "; ".join(categories)
            row["Matched_annotation_terms"] = "; ".join(sorted(set(matched_terms)))
            hits.append(row)

            for category in categories:
                counts[category] += 1

    hits.sort(
        key=lambda r: (
            r["Functional_category_candidate"],
            -float(r["DL_expression_percentile"])
        )
    )

    fields = [
        "DL_ORF", "DL_Average_TPM", "DL_expression_percentile",
        "Functional_category_candidate", "Matched_annotation_terms",
        "KOfam_KO", "KOfam_definition",
        "EggNOG_preferred_name", "EggNOG_description",
        "Direct_Pfam",
        "NI_RBH", "SR_RBH", "PT_RBH", "TP_RBH",
        "RBH_pattern", "RBH_count",
        "NI_protein", "SR_protein", "PT_protein", "TP_protein",
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(hits)

    print("Candidates found:", len(hits))
    for category in PATTERNS:
        print(category, counts[category])
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
