#!/usr/bin/env python3
"""Script 47. Module level dumbbell (main figure): every function in fig_alkaline_modules.png, compared with
MMETSP cultures, plus DUF families with a transcript ranking above cultures and the
housekeeping reference.

One row per module. The row shows the module's most expressed robust transcript
(>= MIN_LIBS libraries, >= MIN_PID median identity), whether or not it ranks above
cultures, so rows are not chosen by the outcome. The filled circle is black when that
transcript was called HIGHER on both metrics and grey otherwise. Two columns on the
right give the libraries for that transcript and, for the whole module, the number of
transcripts (contigs) called higher out of those scored. Modules with no robust
transcript show only the Deer Lake point of their most expressed transcript, in a
separate colour.

Fermentation rows follow fig_alkaline_modules.png and are built from Pfam names;
pyruvate, phosphate dikinase is placed under inorganic carbon, as in that figure.
DUF families are shown only if they contain a transcript above cultures, and their row
shows that transcript.
Bacterial modules removed from the query set (NhaA/NhaB/NhaC, Mrp, microbial
rhodopsin) are not shown.

Inputs: expanded/comparison_per_gene.tsv, alkaline_adaptation_candidates.csv.
Outputs (in expanded/): Fig_MMETSP_dumbbell_modules.png (1000 dpi), dumbbell_modules_rows.tsv.
Lightweight: run interactively.
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt

# ---------------- edit these ----------------
WORK = Path("/work/ebg_lab/eb/diatom_consortia/mmetsp_comparison")   # laptop: folder with the two inputs
DIR = WORK / "expanded"      # folder holding comparison_per_gene.tsv (laptop: may be WORK itself)
ALKALINE_CSV = WORK / "alkaline_adaptation_candidates.csv"
MIN_LIBS, MIN_PID = 20, 50
AXIS_MIN = 0.001
SHOW_DUF_FAMILIES = True     # DUF families with at least one transcript above cultures
NO_ORTHO = "#CC79A7"         # Deer Lake point for functions not compared with cultures
MARK_ABOVE = False           # True adds an asterisk after functions with a transcript above cultures
# fermentation rows: Pfam name -> (row label, section)
FERM = {"ADH_Fe_C": ("Aldehyde/alcohol dehydrogenase", "Fermentation"),
        "Aldedh": ("Aldehyde/alcohol dehydrogenase", "Fermentation"),
        "Fe-ADH": ("Aldehyde/alcohol dehydrogenase", "Fermentation"),
        "PFL-like": ("Pyruvate formate lyase", "Fermentation"),
        "Gly_radical": ("Pyruvate formate lyase", "Fermentation"),
        "Fer4_12": ("Pyruvate formate lyase", "Fermentation"),
        "PEPCK_PPi_lobe_2": ("Pyrophosphate dependent PEP carboxykinase", "Fermentation"),
        "Formate transporter": ("Formate transporter", "Fermentation"),
        "PPDK_N": ("Pyruvate phosphate dikinase", "Inorganic carbon"),
        "PEP-utilizers": ("Pyruvate phosphate dikinase", "Inorganic carbon"),
        "PEP-utilizers_C": ("Pyruvate phosphate dikinase", "Inorganic carbon")}
HK = {"K03231": "Elongation factor 1α", "K03234": "Elongation factor 2", "K05692": "Actin",
      "K07374": "α tubulin", "K07375": "β tubulin", "K08770": "Ubiquitin"}
SECTIONS = {"Inorganic carbon": ("#2A9FD6", "#E3F2FA"),
            "Na⁺/H⁺ and cation/proton antiport": ("#D62728", "#F9E3E3"),
            "Primary pumps and Na⁺ energetics": ("#EE8A12", "#FCF1E1"),
            "Nutrients made scarce by high pH": ("#2E8B57", "#E6F0EA"),
            "Osmolytes": ("#8E44CC", "#F0EAF8"),
            "Fermentation": ("#8B4A3C", "#F2E8E6"),
            "Unknown function": ("#6A5A7A", "#EDEAF0"),
            "Housekeeping reference": ("#777777", "#EEEEEE")}
ORDER = list(SECTIONS)
# --------------------------------------------

ion = lambda t: (str(t).replace("Na+", "Na⁺").replace("H+", "H⁺").replace("K+", "K⁺")
                 .replace("Ca2+", "Ca²⁺").replace("HCO3-", "HCO₃⁻").replace("Cl-", "Cl⁻"))
short = lambda m: ion(re.sub(r"\s*\(.*\)$", "", re.sub(r"^\d+\s*", "", str(m))))

R = pd.read_csv(DIR / "comparison_per_gene.tsv", sep="\t")
need = ["orf", "group", "module", "user_pct", "user_TPM", "ref_pct_median", "n_libs", "median_pident", "call"]
q_lo = [c for c in R.columns if "pct" in c and "25" in c]
q_hi = [c for c in R.columns if "pct" in c and "75" in c]
missing = [c for c in need if c not in R.columns]
if missing or len(q_lo) != 1 or len(q_hi) != 1:
    raise SystemExit(f"Column check failed. Missing: {missing}; IQR columns: {q_lo} {q_hi}\nAvailable: {list(R.columns)}")
Q25, Q75 = q_lo[0], q_hi[0]

a = pd.read_csv(ALKALINE_CSV)
a["module"] = a.module.str.strip()
mech = a.drop_duplicates("module").set_index("module").mechanism.map(short)

R["contig"] = R.orf.str.rsplit(".", n=1).str[0]
R["robust"] = (R.n_libs >= MIN_LIBS) & (R.median_pident >= MIN_PID)
R["higher"] = R.call.astype(str).str.startswith("HIGHER")
R["module"] = R.module.astype(str)

# assign every ORF to its display row(s)
rows = []
for _, r in R.iterrows():
    if r.group == "alkaline":
        for m in r.module.split("; "):
            rows.append((r.name, ion(m), mech.get(m, "Alkaline candidates")))
    elif r.group == "fermentation":
        if r.module in FERM:
            lab, sec = FERM[r.module]
            rows.append((r.name, lab, sec))
    elif r.group == "housekeeping":
        rows.append((r.name, HK.get(r.module, r.module), "Housekeeping reference"))
    elif r.group == "DUF" and SHOW_DUF_FAMILIES:
        rows.append((r.name, r.module, "Unknown function"))
M = pd.DataFrame(rows, columns=["idx", "row", "section"]).join(R, on="idx")
# the alkaline module is called "Pyruvate phosphate dikinase" too, so both sources merge into one row

summ = []
for (sec, row), g in M.groupby(["section", "row"]):
    rob = g[g.robust]
    contigs = rob.groupby("contig").higher.any()
    k, n = int(contigs.sum()), len(contigs)
    if sec == "Unknown function" and k == 0:
        continue                                  # only DUF families with a higher transcript
    # DUF families are shown because they contain a higher transcript, so show that one;
    # every other row shows its most expressed robust transcript, whatever its call
    pool = rob[rob.higher] if sec == "Unknown function" else rob
    rep = pool.sort_values("user_TPM", ascending=False).head(1)
    d = {"section": sec, "row": row, "k_higher": k, "n_scored": n, "n_orfs_all": g.orf.nunique()}
    if len(rep):
        r = rep.iloc[0]
        d.update(orf=r.orf, n_libs=int(r.n_libs), user_pct=r.user_pct, ref_pct_median=r.ref_pct_median,
                 q25=r[Q25], q75=r[Q75], rep_higher=bool(r.higher), compared=True)
    else:
        # no transcript with enough culture orthologs: show the most expressed transcript alone
        r = g.sort_values("user_TPM", ascending=False).iloc[0]
        d.update(orf=r.orf, n_libs=int(r.n_libs) if pd.notna(r.n_libs) else 0,
                 user_pct=r.user_pct, compared=False)
    summ.append(d)
S = pd.DataFrame(summ)
S["s_rank"] = S.section.map({s: i for i, s in enumerate(ORDER)}).fillna(len(ORDER))
S["sort_pct"] = S.user_pct.fillna(-1)
S = S.sort_values(["s_rank", "sort_pct"], ascending=[True, False]).reset_index(drop=True)
S.loc[S.row == "RuBisCO", "row"] = "RuBisCO (plastid encoded)"
S.drop(columns=["s_rank", "sort_pct"]).to_csv(DIR / "dumbbell_modules_rows.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(S[["section", "row", "orf", "compared", "n_libs", "user_pct", "ref_pct_median", "rep_higher", "k_higher", "n_scored"]]
      .to_string(index=False))

# ---------------- figure ----------------
tr = lambda s: np.clip(100 - s, AXIS_MIN, 100)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6})
ys, heads, pos, prev = [], [], 0.0, None
for sec in S.section:
    if sec != prev:
        heads.append([sec, pos])
        pos -= 1.0
        prev = sec
    ys.append(pos)
    pos -= 1.0
S["y"] = ys
S["label"] = [f"{r}*" if (MARK_ABOVE and k > 0) else r for r, k in zip(S.row, S.k_higher)]
fig, ax = plt.subplots(figsize=(6.0, 0.9 + 0.19 * (len(S) + len(heads))))
for sec, top_y in heads:
    txt, band = SECTIONS.get(sec, ("#666666", "#F2F2F2"))
    bot_y = S.y[S.section == sec].min()
    ax.axhspan(bot_y - 0.5, top_y + 0.5, color=band, lw=0, zorder=0)
    ax.text(0.01, top_y, sec, transform=ax.get_yaxis_transform(), color=txt, fontsize=7,
            fontweight="bold", fontstyle="italic", va="center")
for _, r in S.iterrows():
    if not r.compared:
        ax.scatter(tr(r.user_pct), r.y, s=18, c=NO_ORTHO, lw=0, zorder=4, clip_on=False)
    else:
        lo, hi, med, dl = tr(r.q25), tr(r.q75), tr(r.ref_pct_median), tr(r.user_pct)
        ax.plot([lo, hi], [r.y, r.y], c="#9E9E9E", lw=3, solid_capstyle="butt", zorder=1)
        ax.plot([med, dl], [r.y, r.y], c="#333333", lw=0.8, zorder=2)
        ax.scatter(med, r.y, s=18, facecolors="white", edgecolors="#333333", lw=0.8, zorder=3, clip_on=False)
        ax.scatter(dl, r.y, s=18, c="#222222", lw=0, zorder=4, clip_on=False)
        ax.text(1.01, r.y, f"{int(r.n_libs)}", transform=ax.get_yaxis_transform(), va="center",
                fontsize=6, color="#555555")
ax.text(1.01, 0.5, "Libraries", transform=ax.get_yaxis_transform(), va="center", fontsize=6,
        color="#555555", fontweight="bold")
ax.set_yticks(S.y, S.label)
ax.set_ylim(pos + 0.5, 1.0)
ax.set_xscale("log")
ticks = [t for t in [100, 10, 1, 0.1, 0.01, 0.001] if t >= AXIS_MIN]
ax.set_xticks(ticks, [f"{100 - t:g}" for t in ticks])        # tick labels as percentiles
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xlim(100, AXIS_MIN)
ax.set_xlabel("Expression percentile within transcriptome")
ax.spines[["top", "right"]].set_visible(False)
h1 = ax.scatter([], [], s=18, facecolors="white", edgecolors="#333333", lw=0.8, label="MMETSP median")
h2 = ax.plot([], [], c="#9E9E9E", lw=3, label="MMETSP interquartile range")[0]
h3 = ax.scatter([], [], s=18, c="#222222", lw=0, label="Deer Lake")
h4 = ax.scatter([], [], s=18, c=NO_ORTHO, lw=0, label="Deer Lake, not compared with cultures")
import matplotlib.transforms as mtransforms
LEGEND_GAP = 30   # points between the bottom of the axis and the top of the legend
ax.legend(handles=[h1, h2, h3, h4], frameon=False, ncol=2, fontsize=6.5, loc="upper center",
          bbox_to_anchor=(0.5, 0),
          bbox_transform=mtransforms.offset_copy(ax.transAxes, fig=fig, y=-LEGEND_GAP, units="points"))
fig.savefig(DIR / "Fig_MMETSP_dumbbell_modules.png", dpi=1000, bbox_inches="tight")
print("Wrote", DIR / "Fig_MMETSP_dumbbell_modules.png", "and dumbbell_modules_rows.tsv")
plt.show()
