#!/usr/bin/env python3
"""Prepared plotting script for DL expression percentile + four-reference RBH matrix."""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", default="DL_expression_RBH_matrix")
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    df["DL_expression_percentile"] = pd.to_numeric(df["DL_expression_percentile"])

    rename = {
        "Photosystem II D1 protein (PsbA)": "PsbA",
        "Photosystem I P700 apoprotein A1 (PsaA)": "PsaA",
        "Fucoxanthin chlorophyll a/c binding protein": "Fucoxanthin Chl a/c protein",
        "Carbonic anhydrase, CA-like family": "Carbonic anhydrase (CA-like)",
        "Sedoheptulose-bisphosphatase": "Sedoheptulose bisphosphatase",
        "V-type H+ ATPase proteolipid subunit": "V-type H+ ATPase",
        "Calcium-activated chloride channel": "Ca2+-activated Cl- channel",
        "SLC9 sodium/hydrogen exchanger": "SLC9 Na+/H+ exchanger",
        "L-ascorbate peroxidase": "Ascorbate peroxidase",
        "Putative silicon transporter family protein": "Putative silicon transporter",
    }

    df["label"] = df["Curated_function"].map(lambda x: rename.get(x, x))

    labels = []
    for _, row in df.iterrows():
        label = row["label"]
        if (df["label"] == label).sum() > 1:
            labels.append(f"{label} ({row['DL_ORF']})")
        else:
            labels.append(label)
    df["plot_label"] = labels

    category_names = {
        "Photosynthesis": "Photosynthesis",
        "Carbon fixation & CCM": "Carbon fixation and CCM",
        "Silica metabolism": "Silica metabolism",
        "Ion homeostasis / osmoregulation": "Ion homeostasis",
        "Urea & nitrogen metabolism": "Nitrogen and urea metabolism",
        "Oxidative stress": "Oxidative stress",
    }

    category_colors = {
        "Photosynthesis": "#0072B2",
        "Carbon fixation & CCM": "#009E73",
        "Silica metabolism": "#56B4E9",
        "Ion homeostasis / osmoregulation": "#E69F00",
        "Urea & nitrogen metabolism": "#CC79A7",
        "Oxidative stress": "#D55E00",
    }

    n = len(df)
    fig = plt.figure(figsize=(9.2, max(10.2, n * 0.38)))

    ax_label = fig.add_axes([0.05, 0.07, 0.43, 0.87])
    ax_expr = fig.add_axes([0.49, 0.07, 0.16, 0.87])
    ax_rbh = fig.add_axes([0.69, 0.07, 0.20, 0.87])

    for ax in (ax_label, ax_expr, ax_rbh):
        ax.set_ylim(n - 0.5, -0.5)

    ax_label.set_xlim(0, 1)
    ax_label.axis("off")

    for i, row in df.iterrows():
        ax_label.add_patch(Rectangle(
            (0.00, i - 0.45), 0.025, 0.90,
            facecolor=category_colors[row["Functional_category"]],
            edgecolor="none"
        ))
        ax_label.text(0.045, i, row["plot_label"],
                      ha="left", va="center", fontsize=8.1)

    norm = Normalize(vmin=90, vmax=100)
    cmap = plt.get_cmap("viridis")

    for i, pctl in enumerate(df["DL_expression_percentile"]):
        ax_expr.add_patch(Rectangle(
            (0, i - 0.39), 1, 0.78,
            facecolor=cmap(norm(pctl)),
            edgecolor="white", linewidth=0.6
        ))
        ax_expr.text(
            0.5, i, f"{pctl:.2f}",
            ha="center", va="center", fontsize=7.8,
            color="white" if pctl >= 96 else "black"
        )

    ax_expr.set_xlim(0, 1)
    ax_expr.set_xticks([0.5])
    ax_expr.set_xticklabels(["DL expression\npercentile"], fontsize=9)
    ax_expr.xaxis.tick_top()
    ax_expr.set_yticks([])
    for spine in ax_expr.spines.values():
        spine.set_visible(False)

    species = [
        ("NI_RBH", "NI"),
        ("SR_RBH", "SR"),
        ("PT_RBH", "PT"),
        ("TP_RBH", "TP"),
    ]

    ax_rbh.set_xlim(-0.5, 3.5)
    ax_rbh.set_xticks(range(4))
    ax_rbh.set_xticklabels([label for _, label in species], fontsize=9)
    ax_rbh.xaxis.tick_top()
    ax_rbh.set_yticks([])

    for i, row in df.iterrows():
        for j, (column, _) in enumerate(species):
            value = str(row[column])

            ax_rbh.add_patch(Rectangle(
                (j - 0.45, i - 0.39), 0.90, 0.78,
                facecolor="white", edgecolor="#D0D0D0", linewidth=0.6
            ))

            if value == "1":
                ax_rbh.scatter(j, i, s=65, facecolor="black",
                               edgecolor="black", linewidth=0.6, zorder=3)
            elif value == "0":
                ax_rbh.scatter(j, i, s=65, facecolor="white",
                               edgecolor="#777777", linewidth=1.0, zorder=3)
            else:
                ax_rbh.add_patch(Rectangle(
                    (j - 0.45, i - 0.39), 0.90, 0.78,
                    facecolor="#E6E6E6", edgecolor="#D0D0D0", linewidth=0.6
                ))
                ax_rbh.text(j, i, "-", ha="center", va="center",
                            fontsize=10, color="#666666")

    for spine in ax_rbh.spines.values():
        spine.set_visible(False)

    for category in dict.fromkeys(df["Functional_category"]):
        idx = df.index[df["Functional_category"] == category].tolist()
        start, end = min(idx), max(idx)

        if end < n - 1:
            y = end + 0.5
            for ax in (ax_label, ax_expr, ax_rbh):
                ax.axhline(y, color="#BDBDBD", linewidth=0.8)

        ax_label.text(
            0.00, start - 0.48,
            category_names[category],
            ha="left", va="bottom",
            fontsize=8.5, fontweight="bold",
            color=category_colors[category]
        )

    ax_rbh.text(
        1.5, -1.22, "Reciprocal best protein hit",
        ha="center", va="bottom",
        fontsize=9.5, fontweight="bold"
    )

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none",
                   markerfacecolor="black", markeredgecolor="black",
                   markersize=6, label="RBH detected"),
        plt.Line2D([0], [0], marker="o", linestyle="none",
                   markerfacecolor="white", markeredgecolor="#777777",
                   markersize=6, label="No RBH detected"),
        Patch(facecolor="#E6E6E6", edgecolor="#D0D0D0",
              label="Not assessed"),
    ]

    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.69, 0.005), ncol=3,
               frameon=False, fontsize=8)

    cax = fig.add_axes([0.49, 0.025, 0.16, 0.012])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_ticks([90, 95, 100])
    cbar.ax.tick_params(labelsize=7, length=2)
    cbar.set_label("DL expression percentile", fontsize=7)

    for ext in ("png", "pdf", "svg"):
        path = f"{args.out_prefix}.{ext}"
        if ext == "png":
            plt.savefig(path, dpi=1000, bbox_inches="tight")
        else:
            plt.savefig(path, bbox_inches="tight")
        print("Saved:", path)

if __name__ == "__main__":
    main()
