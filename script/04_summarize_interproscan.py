#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path.cwd()
protein_fasta = base / "01_input/diatom_predicted_proteins.fa"
interpro_file = base / "05_interproscan/DL_diatom_braker4_ET_interproscan.tsv"
outdir = base / "06_combined_annotation"
outdir.mkdir(parents=True, exist_ok=True)

summary_by_protein_out = outdir / "DL_diatom_interproscan_summary_by_protein.tsv"
all_proteins_out = outdir / "DL_diatom_all_proteins_with_interproscan_summary.tsv"
counts_out = outdir / "DL_diatom_interproscan_analysis_counts.tsv"
stats_out = outdir / "DL_diatom_interproscan_summary_stats.txt"

def read_fasta_ids(path):
    ids = []
    lengths = []
    current = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    ids.append(current)
                    lengths.append(len("".join(seq)))
                current = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if current is not None:
            ids.append(current)
            lengths.append(len("".join(seq)))
    return pd.DataFrame({"protein_id": ids, "protein_length": lengths})

def join_unique(series):
    vals = []
    for x in series.dropna().astype(str):
        x = x.strip()
        if x and x != "-":
            vals.extend([v.strip() for v in x.split("|") if v.strip() and v.strip() != "-"])
    return "; ".join(sorted(set(vals)))

proteins = read_fasta_ids(protein_fasta)
cols = [
    "protein_id", "md5", "length", "analysis", "signature_accession",
    "signature_description", "start", "end", "score", "status", "date",
    "interpro_accession", "interpro_description", "go_terms", "pathway_annotations"
]
ips = pd.read_csv(interpro_file, sep="\t", names=cols, dtype=str, keep_default_na=False)

counts = ips["analysis"].value_counts().reset_index()
counts.columns = ["analysis", "raw_rows"]
counts.to_csv(counts_out, sep="\t", index=False)

summary = ips.groupby("protein_id", as_index=False).agg(
    interproscan_raw_rows=("protein_id", "size"),
    interproscan_analyses=("analysis", join_unique),
    signature_accessions=("signature_accession", join_unique),
    signature_descriptions=("signature_description", join_unique),
    interpro_accessions=("interpro_accession", join_unique),
    interpro_descriptions=("interpro_description", join_unique),
    go_terms=("go_terms", join_unique),
    pathway_annotations=("pathway_annotations", join_unique),
)
summary.to_csv(summary_by_protein_out, sep="\t", index=False)

all_proteins = proteins.merge(summary, on="protein_id", how="left")
all_proteins["has_interproscan_hit"] = all_proteins["interproscan_raw_rows"].notna().map({True: "yes", False: "no"})
all_proteins.to_csv(all_proteins_out, sep="\t", index=False)

with open(stats_out, "w") as out:
    out.write(f"Raw InterProScan rows: {len(ips)}\n")
    out.write(f"Proteins with InterProScan hits: {summary['protein_id'].nunique()} / {len(proteins)}\n")
    out.write(f"Percent with InterProScan hits: {summary['protein_id'].nunique() / len(proteins) * 100:.2f}%\n")

print("Wrote:", summary_by_protein_out)
print("Wrote:", all_proteins_out)
print("Wrote:", counts_out)
print("Wrote:", stats_out)
