#!/usr/bin/env python3
"""Create BLASTN-unmatched Deer Lake gene subsets without a length cutoff."""

import argparse
from pathlib import Path
import pandas as pd


def read_fasta(path):
    records = {}
    h = None
    seq = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if h is not None:
                    records[h] = ''.join(seq)
                h = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if h is not None:
            records[h] = ''.join(seq)
    return records


def informative_annotation(x):
    s = '' if pd.isna(x) else str(x).strip().lower()
    if not s or s == 'unannotated':
        return False
    return ('uncharacterized' not in s) and ('hypothetical' not in s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--master', required=True)
    p.add_argument('--query-fasta', required=True)
    p.add_argument('--out-prefix', default='DL')
    args = p.parse_args()

    df = pd.read_csv(args.master, sep='\t', dtype=str)
    hit_cols = ['NI_hit','SR_hit','PT_hit','TP_hit']
    for c in hit_cols:
        if c not in df.columns:
            raise SystemExit(f'ERROR: missing {c}')

    nohit = df[df[hit_cols].fillna('no').eq('no').all(axis=1)].copy()
    tpm = pd.to_numeric(nohit['diatom_Average_TPM'], errors='coerce')
    expressed = nohit[tpm > 0].copy()
    expressed['_TPM_numeric'] = pd.to_numeric(expressed['diatom_Average_TPM'], errors='coerce')
    expressed = expressed.sort_values('_TPM_numeric', ascending=False).drop(columns='_TPM_numeric')
    annotated = expressed[expressed['functional_annotation'].map(informative_annotation)].copy()

    prefix = args.out_prefix
    nohit.to_csv(f'{prefix}_no_comparator_hit.tsv', sep='\t', index=False, na_rep='')
    Path(f'{prefix}_no_comparator_hit_gene_roots.txt').write_text('\n'.join(nohit['gene_root'].astype(str)) + '\n')
    expressed.to_csv(f'{prefix}_no_comparator_hit_expressed_sorted_by_TPM.tsv', sep='\t', index=False, na_rep='')
    annotated.to_csv(f'{prefix}_unmatched_expressed_annotated_sorted_TPM.tsv', sep='\t', index=False, na_rep='')

    seqs = read_fasta(args.query_fasta)
    with open(f'{prefix}_no_comparator_hit_genes.fasta', 'w') as out:
        for root in nohit['gene_root'].astype(str):
            if root not in seqs:
                raise SystemExit(f'ERROR: {root} missing from query FASTA')
            seq = seqs[root]
            out.write(f'>{root}\n')
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + '\n')

    print('No comparator hit:', len(nohit))
    print('No comparator hit + TPM > 0:', len(expressed))
    print('No comparator hit + TPM > 0 + informative annotation:', len(annotated))


if __name__ == '__main__':
    main()
