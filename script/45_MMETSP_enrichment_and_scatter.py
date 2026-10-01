#!/usr/bin/env python3
"""Script 45. Enrichment of 'higher than cultures' genes and the all-gene scatter figure.

Reads expanded/comparison_per_gene.tsv written by 03_compare.py.

Scored genes: at least MIN_LIBS libraries and median identity >= MIN_PID (the robust set).
A gene is 'higher' when 03_compare.py called it HIGHER on both metrics.
For every module (alkaline modules, fermentation genes, themes) and every Pfam family
among the DUFs, a one-sided Fisher exact test asks whether higher genes are
overrepresented relative to all other scored non-housekeeping genes; P values are
corrected with Benjamini-Hochberg.

Outputs (in expanded/):
  enrichment_modules.tsv
  Fig_MMETSP_all_genes_scatter.png   (1000 dpi)
Lightweight: run interactively.
"""
from math import lgamma, exp
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- edit these ----------------
DIR = Path("/work/ebg_lab/eb/diatom_consortia/mmetsp_comparison/expanded")
MIN_LIBS = 20
MIN_PID = 50
MIN_GENES = 3          # modules or families with fewer scored genes are not tested
FDR = 0.05
HK_DOMINANT_ONLY = True  # plot only the most expressed ORF per housekeeping KO (fragments hidden)
TECHNICAL = ["RuBisCO"]
AXIS_MIN = 0.01          # rank from top (%) at the upper end of both axes
# labels for higher genes, keyed by module (Pfam name for fermentation ORFs);
# one label per contig, placed on its highest ORF
LABELS = {}             # gene labels are left to the dumbbell figure; add {module: label} to label points
OFFSETS = {}             # nudge overlapping labels, e.g. {"NODE_274": (-25, 4)} in points  # plastid encoded: depleted in poly-A MMETSP libraries, so set aside
# --------------------------------------------

R = pd.read_csv(DIR / "comparison_per_gene.tsv", sep="\t")
R["scored"] = (R.n_libs >= MIN_LIBS) & (R.median_pident >= MIN_PID)
R["higher"] = R.call.astype(str).str.startswith("HIGHER")
R["technical"] = R.module.astype(str).str.split("; ").apply(lambda m: any(t in m for t in TECHNICAL))
R["contig"] = R.orf.str.rsplit(".", n=1).str[0]
S = R[R.scored & (R.group != "housekeeping") & ~R.technical].copy()

# one module per gene; alkaline ORFs in several modules count in each
S["test_unit"] = S.module.astype(str).str.split("; ")
S = S.explode("test_unit")
S["test_unit"] = S.group + ": " + S.test_unit


def lchoose(n, k):
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def fisher_greater(a, b, c, d):
    """One-sided Fisher exact P for enrichment in the 2x2 table [[a, b], [c, d]]."""
    n1, n2, k, n = a + b, c + d, a + c, a + b + c + d
    hi = min(n1, k)
    return min(1.0, sum(exp(lchoose(n1, i) + lchoose(n2, k - i) - lchoose(n, k)) for i in range(a, hi + 1)))


genes = S.drop_duplicates("orf")
N, K = len(genes), int(genes.higher.sum())
rows = []
for unit, g in S.groupby("test_unit"):
    n, x = g.orf.nunique(), int(g.drop_duplicates("orf").higher.sum())
    if n < MIN_GENES:
        continue
    rest_n = N - n
    rest_x = K - x
    rows.append((unit, n, x, round(100 * x / n, 1), fisher_greater(x, n - x, rest_x, rest_n - rest_x)))
E = pd.DataFrame(rows, columns=["module", "n_scored", "n_higher", "pct_higher", "p"])
E = E.sort_values("p").reset_index(drop=True)
m = len(E)
E["q_BH"] = (E.p * m / (np.arange(m) + 1))[::-1].cummin()[::-1].clip(upper=1)
E["enriched"] = E.q_BH < FDR

# group level tests (whole modules such as fermentation), per ORF and per contig;
# ORFs from one contig are not independent, so the contig test is the conservative one
G = S.drop_duplicates("orf")
C = G.groupby(["contig"]).agg(group=("group", "first"), higher=("higher", "any")).reset_index()
grows = []
for grp in sorted(G.group.unique()):
    for level, T in (("ORF", G), ("contig", C)):
        n, x = int((T.group == grp).sum()), int(T.higher[T.group == grp].sum())
        rn, rx = len(T) - n, int(T.higher.sum()) - x
        grows.append((f"group: {grp} [{level}]", n, x, round(100 * x / n, 1),
                      fisher_greater(x, n - x, rx, rn - rx)))
Eg = pd.DataFrame(grows, columns=E.columns[:5])
Eg["q_BH"] = np.nan
Eg["enriched"] = np.nan
E = pd.concat([Eg, E], ignore_index=True)
E.to_csv(DIR / "enrichment_modules.tsv", sep="\t", index=False)
print("Group level (one sided Fisher test against all other scored genes):")
print(Eg.iloc[:, :5].to_string(index=False))
print()
tech = R[R.scored & R.technical]
print(f"Set aside as technical ({', '.join(TECHNICAL)}): {len(tech)} scored, {int(tech.higher.sum())} higher")

print(f"Scored non-housekeeping genes (technical set aside): {N:,}; higher: {K:,} ({100 * K / N:.1f}%)")
print(f"Tested modules/families: {m}; enriched at FDR {FDR}: {int((E.enriched == True).sum())}")
print(E[E.module.str.startswith("group:") == False].head(20).to_string(index=False))

# ---------------- scatter figure ----------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6})
top = lambda pct: np.clip(100 - pct, AXIS_MIN, 100)     # rank from the top, in %
A = R[R.scored].copy()
A["x"], A["y"] = top(A.ref_pct_median), top(A.user_pct)

COL = {"fermentation": "#D55E00", "alkaline": "#0072B2", "DUF": "#009E73", "theme": "#CC79A7"}
fig, ax = plt.subplots(figsize=(3.5, 4.1))
base = A[~A.higher & (A.group != "housekeeping") & ~A.technical]
ax.scatter(base.x, base.y, s=2, c="#BBBBBB", lw=0, alpha=0.5, rasterized=True, clip_on=False, label="Other scored genes")
for grp, c in COL.items():
    g = A[A.higher & (A.group == grp) & ~A.technical]
    if len(g):
        ax.scatter(g.x, g.y, s=5, c=c, lw=0, clip_on=False, zorder=3, label=f"Higher: {grp} (n = {len(g)})")
tp = A[A.technical]
ax.scatter(tp.x, tp.y, s=10, marker="D", facecolors="none", edgecolors="#555555", lw=0.5, clip_on=False,
           label="RuBisCO (plastid encoded)")
hkp = A[A.group == "housekeeping"]
if HK_DOMINANT_ONLY:
    hkp = hkp.sort_values("user_TPM", ascending=False).drop_duplicates("module")
ax.scatter(hkp.x, hkp.y, s=10, facecolors="none", edgecolors="black", lw=0.5, clip_on=False, label="Housekeeping" + (" (most expressed copy)" if HK_DOMINANT_ONLY else ""))
lim = (100, AXIS_MIN)
ax.plot([AXIS_MIN, 100], [AXIS_MIN, 100], c="black", lw=0.5, ls="--", zorder=0)
ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim)
L = A[A.higher & ~A.technical & A.module.isin(LABELS)].sort_values("y").drop_duplicates("contig")
for _, r in L.iterrows():
    ax.annotate(LABELS[r.module], (r.x, r.y), xytext=OFFSETS.get(r.contig, (3, 1.5)), textcoords="offset points",
                fontsize=4.5, va="bottom", clip_on=False)
ticks = [t for t in [100, 10, 1, 0.1, 0.01, 0.001] if t >= AXIS_MIN]
lab = [f"{t:g}" for t in ticks]
ax.set_xticks(ticks, lab)
ax.set_yticks(ticks, lab)
ax.set_xlabel("MMETSP cultures, median rank from top (%)")
ax.set_ylabel("Deer Lake, rank from top (%)")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.17),
          ncol=2, markerscale=1.5, handletextpad=0.3, columnspacing=1.0)
fig.tight_layout()
fig.savefig(DIR / "Fig_MMETSP_all_genes_scatter.png", dpi=1000)
print("Wrote", DIR / "Fig_MMETSP_all_genes_scatter.png", "and enrichment_modules.tsv")
