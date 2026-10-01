#!/usr/bin/env python3
"""Script 46. Transcript level dumbbell figure for the expanded MMETSP comparison, in the style of
Fig_MMETSP_robust_dumbbell.png (coloured sections, library counts on the right).

Rows: every robust 'higher' transcript (one row per contig, its most expressed ORF),
one RuBisCO row (plastid encoded, so shown for reference only), and the most
expressed housekeeping copy per KO as the reference.

Alkaline genes are grouped by the 'mechanism' column of alkaline_adaptation_candidates.csv.
Nothing is recomputed: values come from comparison_per_gene.tsv (03_compare.py).
The script stops, listing the columns, if a column it needs is missing.

Outputs (in DIR): Fig_MMETSP_dumbbell_expanded.png (1000 dpi) and dumbbell_rows.tsv.
Lightweight: run interactively.
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- edit these ----------------
WORK = Path("/work/ebg_lab/eb/diatom_consortia/mmetsp_comparison")
DIR = WORK / "expanded"
ALKALINE_CSV = WORK / "alkaline_adaptation_candidates.csv"
MIN_LIBS, MIN_PID = 20, 50           # same robust filter as enrichment_and_scatter.py
AXIS = "log"                         # "log" = rank from top (%), "percentile" = 0-100 as in the original figure
AXIS_MIN = 0.001                     # right end of the log axis
LABELS = {"PEPCK_PPi_lobe_2": "Pyrophosphate dependent PEP carboxykinase",
          "PPDK_N": "Pyruvate, phosphate dikinase", "PEP-utilizers": "Pyruvate, phosphate dikinase",
          "ADH_Fe_C": "Alcohol dehydrogenase (iron containing)",
          "Fer4_12": "Pyruvate formate lyase",   # PFL-like PF02901 and Gly_radical PF01228 confirmed
          "Formate transporter": "Formate transporter",
          "PEP carboxykinase": "Phosphoenolpyruvate carboxykinase",
          "Proline synthesis (P5CS)": "Δ$^1$ pyrroline 5 carboxylate synthetase",
          "Na+-coupled solute symporters (SSS/SLC5/SLC6/SLC13)": "Sodium coupled solute symporter",
          "RuBisCO": "RuBisCO (plastid encoded)",
          "K03231": "Elongation factor 1α", "K03234": "Elongation factor 2", "K05692": "Actin",
          "K07374": "α tubulin", "K07375": "β tubulin", "K08770": "Ubiquitin"}
# section colours: (text colour, band colour); unknown alkaline mechanisms get grey
SECTIONS = {"Inorganic carbon": ("#2A9FD6", "#E3F2FA"),
            "Na⁺/H⁺ and cation/proton antiport": ("#D62728", "#F9E3E3"),
            "Primary pumps and Na⁺ energetics": ("#EE8A12", "#FCF1E1"),
            "Nutrients made scarce by high pH": ("#2E8B57", "#E6F0EA"),
            "Osmolytes": ("#8E44CC", "#F0EAF8"),
            "Fermentation": ("#8B4A3C", "#F2E8E6"),
            "Unknown function": ("#6A5A7A", "#EDEAF0"),
            "Housekeeping reference": ("#777777", "#EEEEEE")}
# --------------------------------------------

R = pd.read_csv(DIR / "comparison_per_gene.tsv", sep="\t")
need = ["orf", "group", "module", "user_pct", "user_TPM", "ref_pct_median", "n_libs", "median_pident", "call"]
q_lo = [c for c in R.columns if "pct" in c and "25" in c]
q_hi = [c for c in R.columns if "pct" in c and "75" in c]
missing = [c for c in need if c not in R.columns]
if missing or len(q_lo) != 1 or len(q_hi) != 1:
    raise SystemExit(f"Column check failed. Missing: {missing}; IQR columns found: {q_lo} {q_hi}\n"
                     f"Available: {list(R.columns)}")
Q25, Q75 = q_lo[0], q_hi[0]
print("Interquartile range columns:", Q25, Q75)

# module -> mechanism, shortened: "1 Inorganic carbon (HCO3- use ...)" -> "Inorganic carbon"
a = pd.read_csv(ALKALINE_CSV)
mech = a.drop_duplicates("module").set_index(a.drop_duplicates("module").module.str.strip()).mechanism
short = lambda m: re.sub(r"\s*\(.*\)$", "", re.sub(r"^\d+\s*", "", str(m))).replace("Na+", "Na⁺").replace("H+", "H⁺")
mech_order = {short(m): int(re.match(r"\d+", str(m)).group()) for m in a.mechanism.dropna().unique()}

R["contig"] = R.orf.str.rsplit(".", n=1).str[0]
R["robust"] = (R.n_libs >= MIN_LIBS) & (R.median_pident >= MIN_PID)
R["higher"] = R.call.astype(str).str.startswith("HIGHER")
R["first_module"] = R.module.astype(str).str.split("; ").str[0]


def section(r):
    if r.group == "fermentation":
        return "Fermentation"
    if r.group == "DUF":
        return "Unknown function"
    if r.group == "housekeeping":
        return "Housekeeping reference"
    return short(mech.get(r.first_module, "Alkaline candidates"))


def contig_label(g):
    for m in g.first_module:
        if m in LABELS:
            return LABELS[m]
    return g.first_module.iloc[0]


H = R[R.robust & R.higher & (R.group != "housekeeping")].copy()
lab = H.groupby("contig").apply(contig_label, include_groups=False).rename("label")
H = H.sort_values("user_TPM", ascending=False).drop_duplicates("contig").merge(lab, on="contig")
# one RuBisCO row: the transcript with the most MMETSP libraries (ties: most expressed)
rub = H.first_module == "RuBisCO"
if rub.any():
    keep = H[rub].sort_values(["n_libs", "user_TPM"], ascending=False).index[0]
    print(f"RuBisCO: {rub.sum()} higher transcripts, showing {H.loc[keep, 'orf']} ({int(H.loc[keep, 'n_libs'])} libraries)")
    H = H[~rub | (H.index == keep)]

K = R[R.robust & (R.group == "housekeeping")].sort_values("user_TPM", ascending=False) \
    .drop_duplicates("first_module").copy()
K["label"] = K.first_module.map(LABELS).fillna(K.first_module)

P = pd.concat([H, K])
P["section"] = P.apply(section, axis=1)
rank = {s: i for i, s in enumerate(sorted(mech_order, key=mech_order.get))}
rank.update({"Alkaline candidates": 50, "Fermentation": 60, "Unknown function": 70, "Housekeeping reference": 80})
P["s_rank"] = P.section.map(rank).fillna(55)
P = P.sort_values(["s_rank", "user_pct"], ascending=[True, False]).reset_index(drop=True)

if AXIS == "log":
    tr = lambda s: np.clip(100 - s, AXIS_MIN, 100)
else:
    tr = lambda s: s
P["dl"], P["med"], P["lo"], P["hi"] = tr(P.user_pct), tr(P.ref_pct_median), tr(P[Q25]), tr(P[Q75])
P[["section", "label", "orf", "contig", "n_libs", "median_pident", "user_pct", "ref_pct_median", Q25, Q75]] \
    .to_csv(DIR / "dumbbell_rows.tsv", sep="\t", index=False)
print(P[["section", "label", "orf", "n_libs", "user_pct", "ref_pct_median"]].to_string(index=False))

# ---------------- figure ----------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6})
ys, heads, pos, prev = [], [], 0.0, None
for sec in P.section:
    if sec != prev:
        heads.append([sec, pos, None])
        pos -= 1.0
        prev = sec
    ys.append(pos)
    pos -= 1.0
P["y"] = ys
for h in heads:
    h[2] = P.y[P.section == h[0]].min()          # lowest row of the section

fig, ax = plt.subplots(figsize=(5.6, 0.9 + 0.2 * (len(P) + len(heads))))
x0, x1 = (100, AXIS_MIN) if AXIS == "log" else (0, 100)
for sec, top_y, bot_y in heads:
    txt, band = SECTIONS.get(sec, ("#666666", "#F2F2F2"))
    ax.axhspan(bot_y - 0.5, top_y + 0.5, color=band, lw=0, zorder=0)
    ax.text(0.01, top_y, sec, transform=ax.get_yaxis_transform(), color=txt, fontsize=7,
            fontweight="bold", fontstyle="italic", va="center")
for _, r in P.iterrows():
    ax.plot([r.lo, r.hi], [r.y, r.y], c="#9E9E9E", lw=3, solid_capstyle="butt", zorder=1)
    ax.plot([r.med, r.dl], [r.y, r.y], c="#333333", lw=0.8, zorder=2)
    ax.scatter(r.med, r.y, s=18, facecolors="white", edgecolors="#333333", lw=0.8, zorder=3, clip_on=False)
    ax.scatter(r.dl, r.y, s=18, c="#222222", lw=0, zorder=4, clip_on=False)
    ax.text(1.01, r.y, f"{int(r.n_libs)}", transform=ax.get_yaxis_transform(), va="center",
            fontsize=6, color="#555555")
ax.text(1.01, 0.5, "Libraries", transform=ax.get_yaxis_transform(), va="center",
        fontsize=6, color="#555555", fontweight="bold")
ax.set_yticks(P.y, P.label)
ax.set_ylim(pos + 0.5, 1.0)
if AXIS == "log":
    ax.set_xscale("log")
    ticks = [t for t in [100, 10, 1, 0.1, 0.01, 0.001] if t >= AXIS_MIN]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("Expression rank from top within transcriptome (%)")
else:
    ax.set_xlabel("Expression percentile within transcriptome")
ax.set_xlim(x0, x1)
ax.spines[["top", "right"]].set_visible(False)
h1 = ax.scatter([], [], s=18, facecolors="white", edgecolors="#333333", lw=0.8, label="MMETSP median")
h2 = ax.plot([], [], c="#9E9E9E", lw=3, label="MMETSP interquartile range")[0]
h3 = ax.scatter([], [], s=18, c="#222222", lw=0, label="Deer Lake")
fig.legend(handles=[h1, h2, h3], frameon=False, ncol=3, fontsize=6.5, loc="upper center",
           bbox_to_anchor=(0.55, 0.035))
fig.savefig(DIR / "Fig_MMETSP_dumbbell_expanded.png", dpi=1000, bbox_inches="tight")
print("Wrote", DIR / "Fig_MMETSP_dumbbell_expanded.png", "and dumbbell_rows.tsv")
