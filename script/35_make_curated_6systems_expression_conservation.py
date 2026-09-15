#!/usr/bin/env python3
"""Create the manually curated six-system, 36-gene review table."""

import argparse
import csv

SELECTED = {
    "NODE_67.p4": ("Photosynthesis", "Photosystem II D1 protein (PsbA)", "Likely plastid encoded"),
    "NODE_63.p1": ("Photosynthesis", "Photosystem I P700 apoprotein A1 (PsaA)", "Likely plastid encoded"),
    "NODE_63.p2": ("Photosynthesis", "Photosystem I P700 apoprotein A2 (PsaB)", "Likely plastid encoded"),
    "NODE_44.p4": ("Photosynthesis", "Cytochrome b6 (PetB)", "Likely plastid encoded"),
    "NODE_30207.p1": ("Photosynthesis", "Fucoxanthin chlorophyll a/c binding protein", "Nuclear encoded candidate"),
    "NODE_7928.p2": ("Photosynthesis", "Fucoxanthin chlorophyll a/c binding protein", "Nuclear encoded candidate"),
    "NODE_6033.p2": ("Carbon fixation & CCM", "RuBisCO small subunit", "Carbon fixation"),
    "NODE_6033.p1": ("Carbon fixation & CCM", "RuBisCO large subunit", "Likely plastid encoded"),
    "NODE_8468.p1": ("Carbon fixation & CCM", "Carbonic anhydrase, CA-like family", "CCM candidate"),
    "NODE_10495.p2": ("Carbon fixation & CCM", "Carbonic anhydrase", "CCM candidate"),
    "NODE_20723.p1": ("Carbon fixation & CCM", "Phosphoribulokinase", "Calvin cycle"),
    "NODE_18074.p1": ("Carbon fixation & CCM", "Sedoheptulose-bisphosphatase", "Calvin cycle"),
    "NODE_8343.p1": ("Silica metabolism", "Putative silicon transporter family protein", "EggNOG Silicon transporter + PF03842"),
    "NODE_5705.p1": ("Silica metabolism", "Putative silicon transporter family protein", "EggNOG Silicon transporter + PF03842"),
    "NODE_9414.p1": ("Silica metabolism", "Putative silicon transporter family protein", "EggNOG Silicon transporter + PF03842"),
    "NODE_10595.p1": ("Silica metabolism", "Putative silicon transporter family protein", "EggNOG Silicon transporter + PF03842"),
    "NODE_10596.p1": ("Silica metabolism", "Putative silicon transporter family protein", "EggNOG Silicon transporter + PF03842"),
    "NODE_12637.p1": ("Silica metabolism", "Putative silicon transporter family protein", "EggNOG Silicon transporter + PF03842"),
    "NODE_14611.p3": ("Ion homeostasis / osmoregulation", "V-type H+ ATPase proteolipid subunit", "Proton transport"),
    "NODE_1723.p1": ("Ion homeostasis / osmoregulation", "Bestrophin family chloride channel", "Anion transport"),
    "NODE_5391.p2": ("Ion homeostasis / osmoregulation", "Potassium channel", "K+ transport"),
    "NODE_5837.p1": ("Ion homeostasis / osmoregulation", "Calcium-activated chloride channel", "Ca2+/Cl- homeostasis"),
    "NODE_1090.p2": ("Ion homeostasis / osmoregulation", "SLC9 sodium/hydrogen exchanger", "Na+/H+ exchange; PF00999"),
    "NODE_4370.p1": ("Ion homeostasis / osmoregulation", "Na+/H+ antiporter", "Na+/H+ exchange; PF00999"),
    "NODE_971.p2": ("Urea & nitrogen metabolism", "Nitrate reductase", "Nitrate assimilation"),
    "NODE_8419.p1": ("Urea & nitrogen metabolism", "Nitrate/nitrite transporter", "Inorganic nitrogen uptake"),
    "NODE_12386.p1": ("Urea & nitrogen metabolism", "Ammonium transporter", "Ammonium uptake"),
    "NODE_21007.p2": ("Urea & nitrogen metabolism", "Glutamine synthetase", "Nitrogen assimilation"),
    "NODE_12193.p1": ("Urea & nitrogen metabolism", "Glutamate dehydrogenase", "Nitrogen metabolism"),
    "NODE_10703.p1": ("Urea & nitrogen metabolism", "Argininosuccinate synthase", "Urea cycle associated"),
    "NODE_421.p3": ("Oxidative stress", "Thioredoxin 2", "Redox regulation"),
    "NODE_2665.p1": ("Oxidative stress", "L-ascorbate peroxidase", "Peroxide detoxification"),
    "NODE_10018.p1": ("Oxidative stress", "Monothiol glutaredoxin", "Redox regulation"),
    "NODE_25198.p1": ("Oxidative stress", "Peroxiredoxin 5", "Peroxide detoxification"),
    "NODE_23676.p1": ("Oxidative stress", "Glutathione peroxidase", "Peroxide detoxification"),
    "NODE_26832.p1": ("Oxidative stress", "Fe/Mn superoxide dismutase", "Superoxide detoxification"),
}

CATEGORY_ORDER = {
    "Photosynthesis": 1,
    "Carbon fixation & CCM": 2,
    "Silica metabolism": 3,
    "Ion homeostasis / osmoregulation": 4,
    "Urea & nitrogen metabolism": 5,
    "Oxidative stress": 6,
}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = {}
    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["DL_ORF"] in SELECTED:
                rows[row["DL_ORF"]] = row

    missing = sorted(set(SELECTED) - set(rows))
    if missing:
        raise SystemExit("Missing selected ORFs: " + ", ".join(missing))

    fields = [
        "Functional_category", "Curated_function", "Curator_note",
        "DL_ORF", "DL_Average_TPM", "DL_expression_percentile",
        "KOfam_KO", "KOfam_definition",
        "EggNOG_preferred_name", "EggNOG_description", "Direct_Pfam",
        "NI_RBH", "SR_RBH", "PT_RBH", "TP_RBH",
        "RBH_pattern", "RBH_count",
        "NI_protein", "SR_protein", "PT_protein", "TP_protein",
    ]

    output = []
    for orf, (category, function, note) in SELECTED.items():
        row = rows[orf].copy()
        row["Functional_category"] = category
        row["Curated_function"] = function
        row["Curator_note"] = note
        output.append(row)

    output.sort(
        key=lambda r: (
            CATEGORY_ORDER[r["Functional_category"]],
            -float(r["DL_expression_percentile"])
        )
    )

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(output)

    print("Selected genes:", len(output))
    print("Saved:", args.output)

if __name__ == "__main__":
    main()
