#!/usr/bin/env python3
"""Create the final undirected minimum-support Hi-C edge table."""

import argparse
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--min-support', type=int, default=2)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    df['read_pair_count'] = pd.to_numeric(df['read_pair_count'], errors='coerce')
    df = df[df['read_pair_count'].notna()].copy()

    nonself = df[df['forward_contig_id'] != df['reverse_contig_id']].copy()
    kept = nonself[nonself['read_pair_count'] >= args.min_support].copy()

    kept['node1'] = kept[['forward_contig_id','reverse_contig_id']].min(axis=1)
    kept['node2'] = kept[['forward_contig_id','reverse_contig_id']].max(axis=1)

    collapsed = kept.groupby(['node1','node2'], as_index=False)['read_pair_count'].sum()
    collapsed.to_csv(args.output, sep='\t', index=False)

    nodes = set(collapsed['node1']) | set(collapsed['node2'])
    print('Non-self contig pair rows:', len(nonself))
    print(f'Rows with >={args.min_support} contacts before collapse:', len(kept))
    print('Unique undirected edges:', len(collapsed))
    print('Unique contigs in network:', len(nodes))
    if len(collapsed):
        print('Min contact count:', collapsed['read_pair_count'].min())
        print('Max contact count:', collapsed['read_pair_count'].max())
        print('Median contact count:', collapsed['read_pair_count'].median())


if __name__ == '__main__':
    main()
