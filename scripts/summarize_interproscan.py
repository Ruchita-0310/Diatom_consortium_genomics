from pathlib import Path
import pandas as pd

base = Path("/work/ebg_lab/eb/diatom_consortia/functional_annotation_swissprot")

ipr_file = base / "05_interproscan/DL_diatom_braker4_ET_interproscan.tsv"
protein_fasta = base / "01_input/diatom_predicted_proteins.fa"

outdir = base / "06_combined_annotation"
outdir.mkdir(exist_ok=True)

out_file = outdir / "DL_diatom_interproscan_summary_by_protein.tsv"
all_file = outdir / "DL_diatom_all_proteins_with_interproscan_summary.tsv"
summary_file = outdir / "DL_diatom_interproscan_summary_stats.txt"

cols = [
    "protein_id",
    "md5",
    "protein_length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "end",
    "score",
    "status",
    "date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathway_annotations",
]

df = pd.read_csv(
    ipr_file,
    sep="\t",
    names=cols,
    dtype=str,
    keep_default_na=False,
)

def clean_value(x):
    x = str(x).strip()
    if x in ["", "-", "nan", "None"]:
        return None
    return x

def unique_join(values):
    cleaned = []
    seen = set()
    for v in values:
        v = clean_value(v)
        if v is None:
            continue
        if v not in seen:
            cleaned.append(v)
            seen.add(v)
    return ";".join(cleaned)

def split_join(values, sep="|"):
    items = []
    seen = set()
    for value in values:
        value = clean_value(value)
        if value is None:
            continue
        for part in value.split(sep):
            part = clean_value(part)
            if part is None:
                continue
            if part not in seen:
                items.append(part)
                seen.add(part)
    return ";".join(items)

# Create compact signature labels
df["signature_label"] = (
    df["analysis"].astype(str)
    + ":"
    + df["signature_accession"].astype(str)
    + ":"
    + df["signature_description"].astype(str)
)

# Create compact InterPro labels
df["interpro_label"] = (
    df["interpro_accession"].astype(str)
    + ":"
    + df["interpro_description"].astype(str)
)

# Remove empty InterPro labels
df.loc[df["interpro_accession"].isin(["-", "", "nan"]), "interpro_label"] = ""

summary = (
    df.groupby("protein_id")
    .agg(
        protein_length=("protein_length", "first"),
        n_interproscan_rows=("protein_id", "size"),
        analyses=("analysis", unique_join),
        signature_accessions=("signature_accession", unique_join),
        signature_descriptions=("signature_description", unique_join),
        signatures=("signature_label", unique_join),
        interpro_accessions=("interpro_accession", unique_join),
        interpro_descriptions=("interpro_description", unique_join),
        interpro_entries=("interpro_label", unique_join),
        go_terms=("go_terms", split_join),
        pathway_annotations=("pathway_annotations", split_join),
    )
    .reset_index()
)

summary.to_csv(out_file, sep="\t", index=False)

# Keep all predicted proteins, including those without InterProScan hits
protein_ids = []
with open(protein_fasta) as handle:
    for line in handle:
        if line.startswith(">"):
            protein_ids.append(line[1:].strip().split()[0])

all_proteins = pd.DataFrame({"protein_id": protein_ids})
all_with_ipr = all_proteins.merge(summary, on="protein_id", how="left")
all_with_ipr.to_csv(all_file, sep="\t", index=False)

n_total = len(all_proteins)
n_raw_rows = len(df)
n_with_ipr = summary["protein_id"].nunique()

analysis_counts = (
    df["analysis"]
    .value_counts()
    .rename_axis("analysis")
    .reset_index(name="raw_rows")
)

analysis_counts_file = outdir / "DL_diatom_interproscan_analysis_counts.tsv"
analysis_counts.to_csv(analysis_counts_file, sep="\t", index=False)

stats = f"""InterProScan summary

Total predicted proteins: {n_total}
Raw InterProScan rows: {n_raw_rows}
Proteins with at least one InterProScan hit: {n_with_ipr}
Percent with InterProScan hit: {(n_with_ipr / n_total) * 100:.2f}%

Output files:
{out_file}
{all_file}
{analysis_counts_file}
"""

with open(summary_file, "w") as handle:
    handle.write(stats)

print(stats)
