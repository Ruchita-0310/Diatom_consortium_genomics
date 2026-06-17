#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

base = Path("/work/ebg_lab/eb/diatom_consortia/functional_annotation_swissprot")

diamond_file = base / "02_diamond/DL_diatom_braker4_ET_vs_swissprot.tsv"
protein_fasta = base / "01_input/diatom_predicted_proteins.fa"

outdir = base / "03_best_hits"
outdir.mkdir(exist_ok=True)

cols = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "stitle",
]

df = pd.read_csv(diamond_file, sep="\t", names=cols)

# Calculate alignment coverage
df["qcov_percent"] = (df["length"] / df["qlen"]) * 100
df["scov_percent"] = (df["length"] / df["slen"]) * 100

# Parse UniProt subject ID: sp|ACCESSION|ENTRY
def parse_uniprot_id(x):
    parts = str(x).split("|")
    if len(parts) >= 3:
        return parts[0], parts[1], parts[2]
    return "", x, ""

df[["uniprot_review_status", "swissprot_accession", "swissprot_entry"]] = df["sseqid"].apply(
    lambda x: pd.Series(parse_uniprot_id(x))
)

df["uniprot_review_status"] = df["uniprot_review_status"].replace({
    "sp": "Swiss-Prot_reviewed",
    "tr": "TrEMBL_unreviewed"
})

def extract_protein_name(title):
    title = str(title)
    title = re.sub(r"^(sp|tr)\|[^|]+\|[^\s]+\s+", "", title)
    return re.split(r"\sOS=", title)[0]

def extract_organism(title):
    match = re.search(r"\bOS=(.*?)\sOX=", str(title))
    return match.group(1) if match else ""

def extract_taxid(title):
    match = re.search(r"\bOX=(\d+)", str(title))
    return match.group(1) if match else ""

def extract_gene_name(title):
    match = re.search(r"\bGN=(.*?)\sPE=", str(title))
    return match.group(1) if match else ""

def extract_pe(title):
    match = re.search(r"\bPE=(\d+)", str(title))
    return match.group(1) if match else ""

df["swissprot_protein_name"] = df["stitle"].apply(extract_protein_name)
df["swissprot_organism"] = df["stitle"].apply(extract_organism)
df["swissprot_taxid"] = df["stitle"].apply(extract_taxid)
df["swissprot_gene_name"] = df["stitle"].apply(extract_gene_name)
df["protein_existence_PE"] = df["stitle"].apply(extract_pe)

# Select best hit per query
df_sorted = df.sort_values(
    by=["qseqid", "evalue", "bitscore", "qcov_percent", "pident"],
    ascending=[True, True, False, False, False],
)

best = df_sorted.drop_duplicates(subset="qseqid", keep="first").copy()

def assign_confidence(row):
    if row["evalue"] <= 1e-20 and row["qcov_percent"] >= 70 and row["pident"] >= 40:
        return "high"
    elif row["evalue"] <= 1e-10 and row["qcov_percent"] >= 50 and row["pident"] >= 30:
        return "medium"
    elif row["evalue"] <= 1e-5 and row["qcov_percent"] >= 30:
        return "low"
    else:
        return "weak_domain_or_fragment"

best["swissprot_hit_confidence"] = best.apply(assign_confidence, axis=1)

strict = best[
    (best["evalue"] <= 1e-10)
    & (best["qcov_percent"] >= 50)
    & (best["pident"] >= 30)
].copy()

# Keep all predicted proteins, including proteins without Swiss-Prot hits
protein_ids = []
with open(protein_fasta) as handle:
    for line in handle:
        if line.startswith(">"):
            protein_ids.append(line[1:].strip().split()[0])

all_proteins = pd.DataFrame({"qseqid": protein_ids})
all_with_swissprot = all_proteins.merge(best, on="qseqid", how="left")

df.to_csv(outdir / "DL_diatom_swissprot_all_hits_with_coverage.tsv", sep="\t", index=False)
best.to_csv(outdir / "DL_diatom_swissprot_best_hits.tsv", sep="\t", index=False)
strict.to_csv(outdir / "DL_diatom_swissprot_best_hits_strict.tsv", sep="\t", index=False)
all_with_swissprot.to_csv(outdir / "DL_diatom_all_proteins_with_swissprot_annotation.tsv", sep="\t", index=False)

n_total = len(all_proteins)
n_all_hit_lines = len(df)
n_hit_proteins = best["qseqid"].nunique()
n_strict = strict["qseqid"].nunique()

confidence_counts = best["swissprot_hit_confidence"].value_counts().to_string()

summary = f"""Swiss-Prot annotation summary

Total predicted proteins: {n_total}
Total DIAMOND hit lines: {n_all_hit_lines}
Proteins with at least one Swiss-Prot hit: {n_hit_proteins}
Proteins with strict Swiss-Prot hit: {n_strict}

Percent with at least one Swiss-Prot hit: {(n_hit_proteins / n_total) * 100:.2f}%
Percent with strict Swiss-Prot hit: {(n_strict / n_total) * 100:.2f}%

Strict filter:
e-value <= 1e-10
query coverage >= 50%
percent identity >= 30%

Confidence counts:
{confidence_counts}
"""

with open(outdir / "DL_diatom_swissprot_annotation_summary.txt", "w") as handle:
    handle.write(summary)

print(summary)
