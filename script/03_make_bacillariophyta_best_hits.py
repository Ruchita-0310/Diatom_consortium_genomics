#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

base = Path.cwd()
protein_fasta = base / "01_input/diatom_predicted_proteins.fa"
raw_file = base / "02_diamond/DL_diatom_braker4_ET_vs_uniprot_bacillariophyta.tsv"
outdir = base / "03_best_hits"
outdir.mkdir(parents=True, exist_ok=True)

prefix = "bacillariophyta"

all_hits_out = outdir / f"DL_diatom_bacillariophyta_all_hits_with_coverage.tsv"
best_hits_out = outdir / f"DL_diatom_bacillariophyta_best_hits.tsv"
best_hits_strict_out = outdir / f"DL_diatom_bacillariophyta_best_hits_strict.tsv"
all_proteins_out = outdir / f"DL_diatom_all_proteins_with_bacillariophyta_annotation.tsv"
summary_out = outdir / f"DL_diatom_bacillariophyta_annotation_summary.txt"

cols = [
    "protein_id", "subject_id", "pident", "alignment_length", "query_length",
    "subject_length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"
]

def read_fasta_ids(path):
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return pd.DataFrame({"protein_id": ids})

def confidence(row):
    e = row["evalue"]
    qcov = row["qcov_percent"]
    pid = row["pident"]
    if e <= 1e-20 and qcov >= 70 and pid >= 40:
        return "high"
    if e <= 1e-10 and qcov >= 50 and pid >= 30:
        return "medium"
    if e <= 1e-5 and qcov >= 30:
        return "low"
    return "weak_domain_or_fragment"

def parse_uniprot_title(subject_id, stitle):
    stitle = str(stitle)
    subject_id = str(subject_id)
    accession = subject_id.split("|")[1] if "|" in subject_id and len(subject_id.split("|")) > 1 else subject_id
    entry = subject_id.split("|")[2] if "|" in subject_id and len(subject_id.split("|")) > 2 else ""
    protein_name = stitle.split(" OS=")[0].strip()
    protein_name = re.sub(r"^sp\|[^|]+\|[^ ]+\s+", "", protein_name)
    protein_name = re.sub(r"^tr\|[^|]+\|[^ ]+\s+", "", protein_name)
    organism = ""
    gene_name = ""
    taxon_id = ""
    protein_evidence = ""
    m = re.search(r" OS=(.*?) OX=", stitle)
    if m:
        organism = m.group(1)
    m = re.search(r" GN=([^ ]+)", stitle)
    if m:
        gene_name = m.group(1)
    m = re.search(r" OX=([^ ]+)", stitle)
    if m:
        taxon_id = m.group(1)
    m = re.search(r" PE=([^ ]+)", stitle)
    if m:
        protein_evidence = m.group(1)
    return pd.Series({
        f"bacillariophyta_accession": accession,
        f"bacillariophyta_entry": entry,
        f"bacillariophyta_protein_name": protein_name,
        f"bacillariophyta_organism": organism,
        f"bacillariophyta_gene_name": gene_name,
        f"bacillariophyta_taxon_id": taxon_id,
        f"bacillariophyta_protein_evidence": protein_evidence,
    })

proteins = read_fasta_ids(protein_fasta)
raw = pd.read_csv(raw_file, sep="\t", names=cols, low_memory=False)

for c in ["pident", "alignment_length", "query_length", "subject_length", "evalue", "bitscore"]:
    raw[c] = pd.to_numeric(raw[c], errors="coerce")

raw["qcov_percent"] = raw["alignment_length"] / raw["query_length"] * 100
raw["scov_percent"] = raw["alignment_length"] / raw["subject_length"] * 100
raw[f"bacillariophyta_hit_confidence"] = raw.apply(confidence, axis=1)
parsed = raw.apply(lambda r: parse_uniprot_title(r["subject_id"], r["stitle"]), axis=1)
raw = pd.concat([raw, parsed], axis=1)

raw = raw.sort_values(
    by=["protein_id", "evalue", "bitscore", "qcov_percent", "pident"],
    ascending=[True, True, False, False, False],
)
raw.to_csv(all_hits_out, sep="\t", index=False)

best = raw.groupby("protein_id", as_index=False).head(1).copy()
best.to_csv(best_hits_out, sep="\t", index=False)

strict = best[(best["evalue"] <= 1e-10) & (best["qcov_percent"] >= 50) & (best["pident"] >= 30)].copy()
strict.to_csv(best_hits_strict_out, sep="\t", index=False)

keep = [
    "protein_id",
    f"bacillariophyta_accession", f"bacillariophyta_entry", f"bacillariophyta_protein_name",
    f"bacillariophyta_organism", f"bacillariophyta_gene_name",
    "evalue", "bitscore", "pident", "qcov_percent", "scov_percent",
    f"bacillariophyta_hit_confidence"
]
annot = best[keep].copy()
annot = annot.rename(columns={
    "evalue": f"bacillariophyta_evalue",
    "bitscore": f"bacillariophyta_bitscore",
    "pident": f"bacillariophyta_pident",
    "qcov_percent": f"bacillariophyta_qcov_percent",
    "scov_percent": f"bacillariophyta_scov_percent",
})

all_proteins = proteins.merge(annot, on="protein_id", how="left")
all_proteins[f"has_bacillariophyta_hit"] = all_proteins[f"bacillariophyta_accession"].notna().map({True: "yes", False: "no"})
all_proteins.to_csv(all_proteins_out, sep="\t", index=False)

with open(summary_out, "w") as out:
    out.write(f"Total predicted proteins: {len(proteins)}\n")
    out.write(f"Total DIAMOND hit lines: {len(raw)}\n")
    out.write(f"Proteins with at least one bacillariophyta hit: {best['protein_id'].nunique()}\n")
    out.write(f"Proteins with strict bacillariophyta hit: {strict['protein_id'].nunique()}\n")
    out.write("\nConfidence counts:\n")
    out.write(best[f"bacillariophyta_hit_confidence"].value_counts(dropna=False).to_string())
    out.write("\n")

print("Wrote:", all_hits_out)
print("Wrote:", best_hits_out)
print("Wrote:", best_hits_strict_out)
print("Wrote:", all_proteins_out)
print("Wrote:", summary_out)
