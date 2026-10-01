#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path("/work/ebg_lab/eb/diatom_consortia/seasonal_organelle_analysis")
METADATA = BASE / "merged_sample_metadata.tsv"
MAPPING_DIR = BASE / "Nitzschia_18S_mapping"
COUNT_DIR = MAPPING_DIR / "read_counts"
OUT = MAPPING_DIR / "Nitzschia_18S_seasonal_normalized.tsv"

metadata = pd.read_csv(METADATA, sep="\t")

rows = []

for _, meta in metadata.iterrows():
    sample = meta["Sample"]

    count_file = COUNT_DIR / f"{sample}.tsv"
    cov_file = MAPPING_DIR / f"{sample}_covstats.txt"

    counts = pd.read_csv(count_file, sep="\t")
    cov = pd.read_csv(cov_file, sep="\t")

    total_reads = int(counts.loc[0, "Total_rRNA_reads"])
    row = cov.iloc[0]

    plus_reads = int(row["Plus_reads"])
    minus_reads = int(row["Minus_reads"])
    mapped_reads = plus_reads + minus_reads

    mapped_fraction = mapped_reads / total_reads

    rows.append(
        {
            "Sample": sample,
            "Season": meta["Season"],
            "Sampling_type": meta["Sampling_type"],
            "Total_rRNA_reads": total_reads,
            "Nitzschia_18S_reads": mapped_reads,
            "Nitzschia_18S_fraction": mapped_fraction,
            "Nitzschia_18S_percent": mapped_fraction * 100,
            "Nitzschia_18S_RPM": mapped_fraction * 1_000_000,
            "Reference_coverage_percent": float(row["Covered_percent"]),
            "Covered_bases": int(row["Covered_bases"]),
            "Reference_length": int(row["Length"]),
            "Average_depth": float(row["Avg_fold"]),
        }
    )

df = pd.DataFrame(rows)

season_order = {
    "Summer": 1,
    "Fall": 2,
    "Winter": 3,
}

type_order = {
    "Mat": 1,
    "Sediment": 2,
}

df["_season_order"] = df["Season"].map(season_order)
df["_type_order"] = df["Sampling_type"].map(type_order)

df = (
    df.sort_values(["_season_order", "_type_order", "Sample"])
    .drop(columns=["_season_order", "_type_order"])
)

df.to_csv(OUT, sep="\t", index=False)

print(df.to_string(index=False))
print()
print(f"Number of samples: {len(df)}")
print(f"Saved: {OUT}")
