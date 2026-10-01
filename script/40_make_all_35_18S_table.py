#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path("/work/ebg_lab/eb/diatom_consortia/seasonal_organelle_analysis")
MAP = BASE / "Nitzschia_18S_mapping"

MERGED_METADATA = BASE / "merged_sample_metadata.tsv"
ADDITIONAL_METADATA = BASE / "additional_sample_metadata.tsv"

OUT = MAP / "Nitzschia_18S_ALL_35_samples.tsv"

rows = []


def process_sample(
    sample,
    season,
    sampling_type,
    dark_months,
    replicate,
    cov_file,
    count_file,
):
    counts = pd.read_csv(count_file, sep="\t")
    total_reads = int(counts.loc[0, "Total_rRNA_reads"])

    cov = pd.read_csv(cov_file, sep="\t")
    x = cov.iloc[0]

    plus_reads = int(x["Plus_reads"])
    minus_reads = int(x["Minus_reads"])
    mapped_reads = plus_reads + minus_reads

    fraction = mapped_reads / total_reads

    rows.append(
        {
            "Sample": sample,
            "Season": season,
            "Sampling_type": sampling_type,
            "Dark_incubation_months": dark_months,
            "Replicate": replicate,
            "Total_rRNA_reads": total_reads,
            "Nitzschia_18S_reads": mapped_reads,
            "Nitzschia_18S_fraction": fraction,
            "Nitzschia_18S_percent": fraction * 100,
            "Nitzschia_18S_RPM": fraction * 1_000_000,
            "Reference_coverage_percent": float(x["Covered_percent"]),
            "Covered_bases": int(x["Covered_bases"]),
            "Reference_length": int(x["Length"]),
            "Average_depth": float(x["Avg_fold"]),
        }
    )


# Initial 14 merged libraries.
merged = pd.read_csv(MERGED_METADATA, sep="\t")

for _, r in merged.iterrows():
    sample = r["Sample"]

    dark_months = pd.NA if r["Sampling_type"] == "Mat" else 0

    process_sample(
        sample=sample,
        season=r["Season"],
        sampling_type=r["Sampling_type"],
        dark_months=dark_months,
        replicate=pd.NA,
        cov_file=MAP / f"{sample}_covstats.txt",
        count_file=MAP / "read_counts" / f"{sample}.tsv",
    )


# Additional 21 libraries.
additional = pd.read_csv(ADDITIONAL_METADATA, sep="\t")

for _, r in additional.iterrows():
    sample = r["Sample"]

    process_sample(
        sample=sample,
        season=r["Season"],
        sampling_type=r["Sampling_type"],
        dark_months=r["Dark_incubation_months"],
        replicate=r["Replicate"],
        cov_file=MAP / "additional" / "covstats" / f"{sample}_covstats.txt",
        count_file=MAP / "additional" / "read_counts" / f"{sample}.tsv",
    )


df = pd.DataFrame(rows)

season_order = {
    "Spring": 1,
    "Summer": 2,
    "Fall": 3,
    "Winter": 4,
}

type_order = {
    "Mat": 1,
    "Sediment": 2,
}

df["_season"] = df["Season"].map(season_order)
df["_type"] = df["Sampling_type"].map(type_order)

df = (
    df.sort_values(
        [
            "_season",
            "_type",
            "Dark_incubation_months",
            "Replicate",
            "Sample",
        ],
        na_position="first",
    )
    .drop(columns=["_season", "_type"])
)

df.to_csv(OUT, sep="\t", index=False, na_rep="NA")

print(df.to_string(index=False))
print()
print(f"Number of samples: {len(df)}")
print(f"Saved: {OUT}")
