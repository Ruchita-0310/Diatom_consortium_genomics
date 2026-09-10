#!/usr/bin/env python3
import argparse, csv, re
from collections import Counter
from pathlib import Path

STRONG_TE = re.compile(
    r'reverse transcriptase|rna-dependent dna polymerase|transposase|transposon|retrovirus-related|'
    r'retroelement|retrotranspos|integrase|\bgag\b|copia protein|dde tnp4|dde superfamily endonuclease', re.I)
LOW_INFO = re.compile(r'^$|^unannotated$|^uncharacterized protein$|^coil$|duf6818', re.I)
CLASSIFIED_TE_PREFIXES = ('LTR/','LINE/','DNA/')


def main():
    p=argparse.ArgumentParser(description='Conservatively classify Deer Lake repeat-overlapping protein models.')
    p.add_argument('--candidates', required=True)
    p.add_argument('--outdir', required=True)
    args=p.parse_args()
    outdir=Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    with open(args.candidates) as fh:
        rows=list(csv.DictReader(fh,delimiter='\t'))
    for r in rows:
        pct=float(r['interspersed_repeat_overlap_percent'])
        ann=(r.get('functional_annotation') or '').strip()
        rep=(r.get('dominant_repeat_class') or '').strip()
        if pct >= 50 and STRONG_TE.search(ann):
            tier='TIER1_high_confidence_TE'
        elif pct >= 75 and rep == 'DNA/MULE-MuDR' and 'SWIM' in ann:
            tier='TIER2_probable_TE_review'
        elif pct >= 75 and rep.startswith(CLASSIFIED_TE_PREFIXES) and LOW_INFO.search(ann):
            tier='TIER2_probable_TE_review'
        else:
            tier='TIER3_retain_review'
        r['TE_filter_tier']=tier

    fields=list(rows[0].keys())
    with open(outdir/'DL_TE_candidate_classification.tsv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)
    for tier in ['TIER1_high_confidence_TE','TIER2_probable_TE_review','TIER3_retain_review']:
        subset=[r for r in rows if r['TE_filter_tier']==tier]
        with open(outdir/f'{tier}.tsv','w',newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(subset)
    tier1=[r for r in rows if r['TE_filter_tier']=='TIER1_high_confidence_TE']
    with open(outdir/'DL_TIER1_TE_gene_roots.txt','w') as fh:
        for r in tier1: fh.write(r['gene_root']+'\n')
    with open(outdir/'DL_TIER1_TE_transcript_ids.txt','w') as fh:
        for r in tier1: fh.write(r['transcript_id']+'\n')

    c=Counter(r['TE_filter_tier'] for r in rows)
    print('Total candidates:', len(rows))
    for tier in ['TIER1_high_confidence_TE','TIER2_probable_TE_review','TIER3_retain_review']:
        print(tier, c[tier])

if __name__=='__main__':
    main()
