#!/usr/bin/env python3
"""Plot Panel A: Deer Lake BLASTN sharing groups versus Average_TPM."""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

COLORS = ['#BDBDBD', '#0072B2', '#009E73', '#E69F00', '#CC79A7']
DOT_COLOR = '#111111'


def clean_yes_no(s):
    return s.astype(str).str.strip().str.lower()


def format_p(p):
    if p == 0:
        return 'p < 1 × 10$^{-300}$'
    if p < 0.001:
        e = int(np.floor(np.log10(p)))
        c = p / (10 ** e)
        return f'p = {c:.2f} × 10$^{{{e}}}$'
    return f'p = {p:.3g}'


def bracket(ax, x1, x2, y, h, text):
    ax.plot([x1,x1,x2,x2], [y,y+h,y+h,y], color='black', lw=1.1)
    ax.text((x1+x2)/2, y+h+0.03, text, ha='center', va='bottom', fontsize=9)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--out-prefix', default='PanelA_BLASTN_TPM')
    args = p.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    for c in ['NI_hit','SR_hit','PT_hit','TP_hit']:
        df[c] = clean_yes_no(df[c])
    df['TPM'] = pd.to_numeric(df['diatom_Average_TPM'], errors='coerce')
    df = df[df['TPM'].notna()].copy()
    df['log10_TPM_plus1'] = np.log10(df['TPM'] + 1)

    nohit = df[(df.NI_hit=='no') & (df.SR_hit=='no') & (df.PT_hit=='no') & (df.TP_hit=='no')]['log10_TPM_plus1']
    ni = df[df.NI_hit=='yes']['log10_TPM_plus1']
    sr = df[df.SR_hit=='yes']['log10_TPM_plus1']
    pt = df[df.PT_hit=='yes']['log10_TPM_plus1']
    tp = df[df.TP_hit=='yes']['log10_TPM_plus1']
    groups = [nohit, ni, sr, pt, tp]
    counts = [len(x) for x in groups]

    labels = [
        f'Deer Lake only\n(n = {counts[0]:,})',
        f'Shared with\n$N.\\ inconspicua$\n(n = {counts[1]:,})',
        f'Shared with\n$S.\\ robusta$\n(n = {counts[2]:,})',
        f'Shared with\n$P.\\ tricornutum$\n(n = {counts[3]:,})',
        f'Shared with\n$T.\\ pseudonana$\n(n = {counts[4]:,})',
    ]

    pvals = [mannwhitneyu(nohit, g, alternative='two-sided').pvalue for g in groups[1:]]
    fig, ax = plt.subplots(figsize=(13,7))
    pos = np.arange(1,6)

    vp = ax.violinplot(groups, positions=pos, widths=0.82, showmeans=False, showmedians=False, showextrema=False, bw_method=0.25)
    for body, color in zip(vp['bodies'], COLORS):
        body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.35); body.set_linewidth(1.0)

    rng = np.random.default_rng(42)
    for x, values in zip(pos, groups):
        jitter = rng.normal(x, 0.065, len(values))
        ax.scatter(jitter, values, s=8, c=DOT_COLOR, alpha=0.35, edgecolors='none', zorder=2, rasterized=True)

    ax.boxplot(groups, positions=pos, widths=0.25, patch_artist=True, showfliers=False,
               medianprops={'color':'white','linewidth':2.0},
               boxprops={'facecolor':'#6E6E6E','edgecolor':'black','linewidth':1.2,'alpha':0.85},
               whiskerprops={'color':'black','linewidth':1.1},
               capprops={'color':'black','linewidth':1.1})

    ymax = max(float(x.max()) for x in groups)
    base_y, step, h = ymax + 0.15, 0.34, 0.08
    for i, pv in enumerate(pvals, start=1):
        bracket(ax, 1, i+1, base_y + (i-1)*step, h, format_p(pv))

    ax.set_xticks(pos); ax.set_xticklabels(labels, fontsize=10.5)
    ax.set_ylabel(r'$\log_{10}(\mathrm{Average\ TPM} + 1)$', fontsize=14)
    ax.tick_params(axis='x', length=0); ax.tick_params(axis='y', labelsize=11)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.text(-0.07, 1.02, 'A', transform=ax.transAxes, fontsize=22, fontweight='bold', ha='left', va='bottom')
    ax.set_ylim(-0.05, base_y + step*len(pvals) + 0.30)
    plt.tight_layout()

    for ext in ['png','pdf','svg']:
        kwargs = {'bbox_inches':'tight'}
        if ext == 'png': kwargs['dpi'] = 600
        plt.savefig(f'{args.out_prefix}.{ext}', **kwargs)
    plt.show()

    print('Group sizes:', counts)
    print('P values vs Deer Lake only:', pvals)


if __name__ == '__main__':
    main()
