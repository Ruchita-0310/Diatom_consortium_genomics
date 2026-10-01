#!/usr/bin/env python3
"""Script 44. Build the expanded query set for the MMETSP comparison.

Input: alkaline_adaptation_candidates.csv, DUF_genes.csv, interesting_genes.csv
and the master expression table.

ORFs are treated as diatom derived. Modules that are bacterial by definition and
themes labelled as non-diatom or viral are dropped (switch off with DROP_NON_DIATOM).
The same 92 housekeeping genes as the 742 gene run are always included, and
percentiles and the housekeeping normalization are computed from the full master
table exactly as before, so results are comparable with the 742 gene run.

Writes OUTDIR/query_genes.tsv, OUTDIR/query_ids.txt and OUTDIR/dropped_orfs.tsv.
Lightweight: run interactively.
"""
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------- edit these ----------------
BASE = Path("/work/ebg_lab/eb/diatom_consortia")
MASTER = BASE / "metatranscriptomics/transdecoder_to_braker_ID_bridge/CLEAN_REBUILD_FROM_RAW/07_expression/master_with_custom_broad_categories.csv"
WORK = BASE / "mmetsp_comparison"
ALKALINE_CSV = WORK / "alkaline_adaptation_candidates.csv"
DUF_CSV = WORK / "DUF_genes.csv"
THEMES_CSV = WORK / "interesting_genes.csv"
OUTDIR = WORK / "expanded"
DROP_NON_DIATOM = True
# ORFs whose group in the CSVs does not match their annotation: orf -> (group, module)
# NODE_2817.p1: KO K21990/K21993 and eggNOG agree on formate transporter (Pfam DUF6113);
# it was in the fermentation group of the 742 gene run
OVERRIDES = {"NODE_2817.p1": ("fermentation", "Formate transporter")}
HIGH_CONFIDENCE_ONLY = True   # alkaline candidates: keep only "high (sources agree)", as in the 742 gene run
# --------------------------------------------

HK = ["K03231", "K03234", "K05692", "K07374", "K07375", "K08770"]  # EF1a, EF2, actin, tubulins, ubiquitin
BACTERIAL_MODULES = ["Bacterial NhaA/NhaB/NhaC",
                     "Mrp multisubunit Na+/H+ antiporter",
                     "Microbial rhodopsin (H+/Na+/Cl- pump or sensory)"]
NON_DIATOM_THEMES = ("A.", "C.", "E.", "G.")   # ssDNA virus, SITs in non-diatoms, dinoflagellate lipase, rhodopsins


def read(path):
    d = pd.read_csv(path, low_memory=False, encoding_errors="replace")
    d.columns = d.columns.str.strip()
    d["orf"] = d["orf"].astype(str).str.strip()
    return d


df = pd.read_csv(MASTER, low_memory=False)
# Average_TPM in the master table is rounded to whole numbers, which creates ties and
# stripes among weakly expressed ORFs. Use the unrounded mean of the replicates instead.
REPS = [c for c in df.columns if c.startswith("TPM_sample")]
if REPS:
    reps = df[REPS].apply(pd.to_numeric, errors="coerce")
    bad = reps.isna().any(axis=1) & df[REPS].notna().any(axis=1)
    if bad.any():
        print(f"WARNING: {bad.sum()} ORFs have non-numeric replicate TPM values; "
              f"they keep the rounded Average_TPM. Examples: {df.orf[bad].head(5).tolist()}")
    exact = reps.mean(axis=1)
    rounded = pd.to_numeric(df.Average_TPM, errors="coerce")
    df["Average_TPM"] = exact.where(~bad & exact.notna(), rounded)
    print("Average_TPM recomputed from", ", ".join(REPS))
else:
    print("WARNING: no TPM_sample columns; using the rounded Average_TPM")
df["pct"] = df.Average_TPM.rank(pct=True) * 100            # rank among all ORFs, as before
hk = df[df.ko.fillna("").str.contains("|".join(HK))]
hk_sum = hk.Average_TPM.sum()

a, d, s = read(ALKALINE_CSV), read(DUF_CSV), read(THEMES_CSV)
s["theme"] = s.theme.astype(str).str.strip()
a["module"] = a.module.astype(str).str.strip()

dropped = []
drop_s = s.theme.str.startswith("A.")                      # viral ORFs are always dropped
drop_d = d["class"].astype(str).str.strip() == "Viral-contig DUF"
drop_a = pd.Series(False, index=a.index)
if HIGH_CONFIDENCE_ONLY:
    drop_a |= ~a.confidence.astype(str).str.startswith("high")
if DROP_NON_DIATOM:
    drop_s |= s.theme.str.startswith(NON_DIATOM_THEMES)
    drop_a |= a.module.isin(BACTERIAL_MODULES)
dropped += [(o, "Sheet3", t) for o, t in zip(s.orf[drop_s], s.theme[drop_s])]
dropped += [(o, "DUF", "Viral-contig DUF") for o in d.orf[drop_d]]
dropped += [(o, "Alkaline", f"{m} [{c}]") for o, m, c in zip(a.orf[drop_a], a.module[drop_a], a.confidence[drop_a])]
s, d, a = s[~drop_s], d[~drop_d], a[~drop_a]
ferm = s.theme.str.startswith("B.")

rows = [(o, "housekeeping", [k for k in HK if k in ko][0], "") for o, ko in zip(hk.orf, hk.ko)]
rows += [(o, "fermentation", p, "") for o, p in zip(s.orf[ferm], s.Pfam_Name[ferm].fillna("unannotated"))]
am = a.groupby("orf").agg(module=("module", lambda v: "; ".join(dict.fromkeys(v))),
                          mech=("mechanism", "first")).reset_index()
rows += [(o, "alkaline", m, str(me).strip()) for o, m, me in zip(am.orf, am.module, am.mech)]
rows += [(o, "theme", t, p if isinstance(p, str) else "")
         for o, t, p in zip(s.orf[~ferm], s.theme[~ferm], s.Pfam_Name[~ferm])]
rows += [(o, "DUF", p, c) for o, p, c in zip(d.orf, d.Pfam_Name, d["class"])]

# first label wins: housekeeping > fermentation > alkaline > theme > DUF
q = pd.DataFrame(rows, columns=["orf", "group", "module", "subclass"]).drop_duplicates("orf")
for o, (g, mod) in OVERRIDES.items():
    hit = q.orf == o
    if hit.any():
        print(f"Override: {o} {q.loc[hit, 'group'].iloc[0]}/{q.loc[hit, 'module'].iloc[0]} -> {g}/{mod}")
        q.loc[hit, ["group", "module", "subclass"]] = [g, mod, "override"]
missing = sorted(set(q.orf) - set(df.orf))
q = q.merge(df[["orf", "Average_TPM", "pct"]], on="orf").rename(columns={"Average_TPM": "user_TPM"})
q["user_pct"] = q.pct.round(3)
q["user_log2_vs_HK"] = np.log2((q.user_TPM + 0.01) / hk_sum)

OUTDIR.mkdir(parents=True, exist_ok=True)
q.drop(columns="pct").to_csv(OUTDIR / "query_genes.tsv", sep="\t", index=False)
q.orf.to_csv(OUTDIR / "query_ids.txt", index=False, header=False)
pd.DataFrame(dropped, columns=["orf", "sheet", "reason"]).to_csv(OUTDIR / "dropped_orfs.tsv", sep="\t", index=False)

print(len(q), "query genes written:", q.group.value_counts().to_dict())
print(len(dropped), "ORF entries dropped (see dropped_orfs.tsv)")
if missing:
    print(f"WARNING: {len(missing)} ORFs not found in the master table, e.g. {missing[:3]}")
