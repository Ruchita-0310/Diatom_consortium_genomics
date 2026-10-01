#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path("/work/ebg_lab/eb/diatom_consortia/seasonal_organelle_analysis")
INPUT = BASE / "Nitzschia_18S_mapping" / "Nitzschia_18S_ALL_35_samples.tsv"
OUT = BASE / "Nitzschia_18S_mapping" / "Nitzschia_18S_group_summary.tsv"

df = pd.read_csv(INPUT, sep="\t", na_values=["NA"])

summary = (
    df.groupby(
        ["Season", "Sampling_type", "Dark_incubation_months"],
        dropna=False,
    )["Nitzschia_18S_percent"]
    .agg(["count", "mean", "std", "min", "max"])
    .reset_index()
)

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

summary["_season"] = summary["Season"].map(season_order)
summary["_type"] = summary["Sampling_type"].map(type_order)

summary = (
    summary.sort_values(
        ["_season", "_type", "Dark_incubation_months"],
        na_position="first",
    )
    .drop(columns=["_season", "_type"])
)

summary.to_csv(OUT, sep="\t", index=False, na_rep="NA")

print(summary.to_string(index=False))
print()
print(f"Saved: {OUT}")
