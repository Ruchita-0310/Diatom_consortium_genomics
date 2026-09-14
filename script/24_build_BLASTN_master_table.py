#!/usr/bin/env python3
"""Build the 14,941-gene Deer Lake BLASTN + TPM master table."""

import argparse
import pandas as pd

BLAST_COLS = [
    'qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend',
    'sstart','send','evalue','bitscore','qlen','slen','qcovs'
]


def fasta_lengths(path):
    lengths = {}
    header = None
    n = 0
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    lengths[header] = n
                header = line[1:].split()[0]
                n = 0
            else:
                n += len(line)
        if header is not None:
            lengths[header] = n
    return lengths


def normalize_root(value):
    return str(value).split('.t', 1)[0]


def best_hits(path, prefix):
    try:
        df = pd.read_csv(path, sep='\t', names=BLAST_COLS, header=None)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=BLAST_COLS)

    if df.empty:
        return pd.DataFrame(columns=[
            'gene_root', f'{prefix}_hit', f'{prefix}_subject', f'{prefix}_pident',
            f'{prefix}_alignment_length', f'{prefix}_qcov', f'{prefix}_evalue', f'{prefix}_bitscore'
        ])

    for c in ['pident','length','evalue','bitscore','qcovs']:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    df['gene_root'] = df['qseqid'].map(normalize_root)
    df = df.sort_values(
        ['gene_root','bitscore','evalue','qcovs','pident','length'],
        ascending=[True,False,True,False,False,False],
        kind='mergesort'
    ).drop_duplicates('gene_root', keep='first')

    out = df[['gene_root','sseqid','pident','length','qcovs','evalue','bitscore']].copy()
    out.insert(1, f'{prefix}_hit', 'yes')
    out = out.rename(columns={
        'sseqid': f'{prefix}_subject',
        'pident': f'{prefix}_pident',
        'length': f'{prefix}_alignment_length',
        'qcovs': f'{prefix}_qcov',
        'evalue': f'{prefix}_evalue',
        'bitscore': f'{prefix}_bitscore',
    })
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--query-fasta', required=True)
    p.add_argument('--metadata-map', required=True)
    p.add_argument('--ni', required=True)
    p.add_argument('--sr', required=True)
    p.add_argument('--pt', required=True)
    p.add_argument('--tp', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--expected', type=int, default=14941)
    args = p.parse_args()

    lengths = fasta_lengths(args.query_fasta)
    master = pd.DataFrame({'gene_root': list(lengths.keys())})
    master['gene_length_bp'] = master['gene_root'].map(lengths)

    meta = pd.read_csv(args.metadata_map, sep='\t', dtype=str)
    if 'gene_root' not in meta.columns:
        if 'gene_id' not in meta.columns:
            raise SystemExit('ERROR: metadata map needs gene_root or gene_id.')
        meta['gene_root'] = meta['gene_id'].map(normalize_root)

    wanted = [
        'gene_root','gene_id','contig_id','diatom_compartment',
        'functional_annotation','diatom_Average_TPM'
    ]
    for col in wanted:
        if col not in meta.columns:
            meta[col] = ''
    meta = meta[wanted].drop_duplicates('gene_root', keep='first')

    master = master.merge(meta, on='gene_root', how='left')
    master = master[[
        'gene_root','gene_id','contig_id','diatom_compartment','gene_length_bp',
        'functional_annotation','diatom_Average_TPM'
    ]]

    for prefix, path in [('NI',args.ni),('SR',args.sr),('PT',args.pt),('TP',args.tp)]:
        hit = best_hits(path, prefix)
        master = master.merge(hit, on='gene_root', how='left')
        master[f'{prefix}_hit'] = master[f'{prefix}_hit'].fillna('no')

    if len(master) != args.expected:
        raise SystemExit(f'ERROR: expected {args.expected} rows, got {len(master)}.')

    master.to_csv(args.out, sep='\t', index=False, na_rep='')

    print('Master rows:', len(master))
    for prefix in ['NI','SR','PT','TP']:
        print(f'{prefix} hit:', int((master[f"{prefix}_hit"] == 'yes').sum()))
    any_hit = master[[f'{x}_hit' for x in ['NI','SR','PT','TP']]].eq('yes').any(axis=1)
    print('Hit in >=1 comparator:', int(any_hit.sum()))
    print('No hit any comparator:', int((~any_hit).sum()))


if __name__ == '__main__':
    main()
