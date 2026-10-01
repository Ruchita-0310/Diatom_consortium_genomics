#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT = (
    "Nitzschia_18S_mapping/"
    "Nitzschia_18S_ALL_35_samples.tsv"
)

OUT_PNG = "Nitzschia_18S_season_dark_incubation.png"

# ============================================================
# READ DATA
# ============================================================

df = pd.read_csv(INPUT, sep="\t", na_values=["NA"])

season_order = ["Spring", "Summer", "Fall", "Winter"]

df["Season"] = pd.Categorical(
    df["Season"],
    categories=season_order,
    ordered=True,
)

# ============================================================
# COLOR BLIND AWARE PALETTE
# ============================================================

season_colors = {
    "Spring": "#009E73",
    "Summer": "#E69F00",
    "Fall": "#D55E00",
    "Winter": "#0072B2",
}

# ============================================================
# STYLE
# ============================================================

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
    }
)

fig, axes = plt.subplots(
    2,
    2,
    figsize=(9.0, 7.0),
)

axA, axB, axC, axD = axes.flatten()

# ============================================================
# PANEL A
# SEASON × SAMPLING CONDITION MATRIX
# ============================================================

conditions = [
    ("Mat", np.nan, "Mat"),
    ("Sediment", 0, "Sediment, 0 months"),
    ("Sediment", 3, "Sediment, 3 months"),
    ("Sediment", 9, "Sediment, 9 months"),
]

condition_labels = [x[2] for x in conditions]

summary_rows = []

for sampling_type, months, label in conditions:
    for season in season_order:
        if sampling_type == "Mat":
            sub = df[
                (df["Sampling_type"] == "Mat")
                & (df["Season"] == season)
            ]
        else:
            sub = df[
                (df["Sampling_type"] == "Sediment")
                & (df["Season"] == season)
                & (df["Dark_incubation_months"] == months)
            ]

        if len(sub) > 0:
            summary_rows.append(
                {
                    "Condition": label,
                    "Season": season,
                    "Mean": sub["Nitzschia_18S_percent"].mean(),
                    "n": len(sub),
                }
            )

summary = pd.DataFrame(summary_rows)

for _, r in summary.iterrows():
    x = season_order.index(r["Season"])
    y = condition_labels.index(r["Condition"])

    # Bubble area is proportional to the mean normalized signal.
    size = 45 + r["Mean"] * 34

    axA.scatter(
        x,
        y,
        s=size,
        color=season_colors[r["Season"]],
        edgecolor="black",
        linewidth=0.6,
        zorder=3,
    )

    axA.text(
        x,
        y,
        f'{r["Mean"]:.1f}',
        ha="center",
        va="center",
        fontsize=7,
        color="black",
        zorder=4,
    )

axA.set_xticks(range(4))
axA.set_xticklabels(season_order)

axA.set_yticks(range(4))
axA.set_yticklabels(condition_labels)

axA.invert_yaxis()

axA.set_xlabel("Season")
axA.set_ylabel("Sampling condition")

axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)

axA.text(
    -0.16,
    1.07,
    "A",
    transform=axA.transAxes,
    fontsize=13,
    fontweight="bold",
)

# ============================================================
# HELPER FOR SEASONAL DOT PLOTS
# ============================================================


def seasonal_dotplot(ax, subset, ylabel):
    for i, season in enumerate(season_order):
        s = subset[subset["Season"] == season].copy()

        if len(s) == 0:
            continue

        values = s["Nitzschia_18S_percent"].values

        if len(values) == 1:
            jitter = np.array([0.0])
        else:
            jitter = np.linspace(-0.055, 0.055, len(values))

        ax.scatter(
            i + jitter,
            values,
            s=52,
            color=season_colors[season],
            edgecolor="black",
            linewidth=0.6,
            zorder=3,
        )

        if len(values) >= 2:
            mean = np.mean(values)
            sd = np.std(values, ddof=1)

            ax.errorbar(
                i,
                mean,
                yerr=sd,
                fmt="_",
                markersize=14,
                markeredgewidth=1.4,
                color="black",
                capsize=3,
                linewidth=1,
                zorder=4,
            )

    ax.set_xticks(range(4))
    ax.set_xticklabels(season_order)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 13)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        direction="out",
        length=3.5,
        width=0.8,
    )

    ax.grid(
        axis="y",
        linewidth=0.4,
        alpha=0.2,
        zorder=0,
    )


# ============================================================
# PANEL B
# MAT ACROSS SEASONS
# ============================================================

mat = df[df["Sampling_type"] == "Mat"].copy()

seasonal_dotplot(
    axB,
    mat,
    r"$\it{Nitzschia}$ like 18S reads (%)",
)

axB.set_xlabel("Season")
axB.set_title("Mat", fontsize=10)

axB.text(
    -0.16,
    1.07,
    "B",
    transform=axB.transAxes,
    fontsize=13,
    fontweight="bold",
)

# ============================================================
# PANEL C
# FRESH SEDIMENT ACROSS SEASONS
# ============================================================

sed0 = df[
    (df["Sampling_type"] == "Sediment")
    & (df["Dark_incubation_months"] == 0)
].copy()

seasonal_dotplot(
    axC,
    sed0,
    r"$\it{Nitzschia}$ like 18S reads (%)",
)

axC.set_xlabel("Season")
axC.set_title("Sediment, 0 months dark incubation", fontsize=10)

axC.text(
    -0.16,
    1.07,
    "C",
    transform=axC.transAxes,
    fontsize=13,
    fontweight="bold",
)

# ============================================================
# PANEL D
# SEDIMENT DARK INCUBATION
# ============================================================

sed = df[
    (df["Sampling_type"] == "Sediment")
    & (df["Season"].isin(["Summer", "Winter"]))
].copy()

months = [0, 3, 9]

offsets = {
    "Summer": -0.10,
    "Winter": 0.10,
}

for season in ["Summer", "Winter"]:
    means = []

    for i, month in enumerate(months):
        s = sed[
            (sed["Season"] == season)
            & (sed["Dark_incubation_months"] == month)
        ]

        values = s["Nitzschia_18S_percent"].values

        if len(values) == 0:
            means.append(np.nan)
            continue

        jitter = np.linspace(
            -0.035,
            0.035,
            len(values),
        )

        xcenter = i + offsets[season]

        axD.scatter(
            xcenter + jitter,
            values,
            s=48,
            color=season_colors[season],
            edgecolor="black",
            linewidth=0.55,
            zorder=4,
        )

        mean = np.mean(values)
        sd = np.std(values, ddof=1)

        means.append(mean)

        axD.errorbar(
            xcenter,
            mean,
            yerr=sd,
            fmt="o",
            markersize=4,
            color="black",
            capsize=3,
            linewidth=1,
            zorder=5,
        )

    axD.plot(
        np.arange(3) + offsets[season],
        means,
        linewidth=1.3,
        color=season_colors[season],
        zorder=2,
    )

axD.set_xticks(range(3))
axD.set_xticklabels(["0", "3", "9"])

axD.set_xlabel("Dark incubation (months)")
axD.set_ylabel(r"$\it{Nitzschia}$ like 18S reads (%)")
axD.set_ylim(0, 13)

axD.spines["top"].set_visible(False)
axD.spines["right"].set_visible(False)

axD.grid(
    axis="y",
    linewidth=0.4,
    alpha=0.2,
    zorder=0,
)

legend = [
    Line2D(
        [0],
        [0],
        marker="o",
        color=season_colors["Summer"],
        label="Summer",
        markersize=5,
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color=season_colors["Winter"],
        label="Winter",
        markersize=5,
    ),
]

axD.legend(
    handles=legend,
    frameon=False,
)

axD.text(
    -0.16,
    1.07,
    "D",
    transform=axD.transAxes,
    fontsize=13,
    fontweight="bold",
)

# ============================================================
# SAVE PNG ONLY
# ============================================================

plt.tight_layout()

plt.savefig(
    OUT_PNG,
    dpi=1000,
    bbox_inches="tight",
)

plt.close()

print(f"Saved: {OUT_PNG}")
