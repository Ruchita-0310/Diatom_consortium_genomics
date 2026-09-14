#!/usr/bin/env python3
"""Prepare the final Deer Lake nucleotide query for the four-genome BLASTN analysis."""

import argparse
from pathlib import Path


def fasta_records(path):
    header = None
    seq = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, ''.join(seq)


def gene_root(identifier):
    # BRAKER representative IDs are typically g123.t1; nucleotide query IDs are g123.
    return identifier.split('.t', 1)[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--gene-fasta', required=True)
    p.add_argument('--final-protein-fasta', required=True)
    p.add_argument('--out-fasta', required=True)
    p.add_argument('--roots-out', required=True)
    p.add_argument('--expected', type=int, default=None)
    args = p.parse_args()

    roots = []
    seen = set()
    for header, _ in fasta_records(args.final_protein_fasta):
        root = gene_root(header)
        if root not in seen:
            roots.append(root)
            seen.add(root)

    Path(args.roots_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.roots_out).write_text('\n'.join(roots) + '\n')

    found = 0
    missing = set(roots)
    with open(args.out_fasta, 'w') as out:
        for header, seq in fasta_records(args.gene_fasta):
            root = gene_root(header)
            if root in seen:
                out.write(f'>{root}\n')
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + '\n')
                found += 1
                missing.discard(root)

    print('Final protein gene roots:', len(roots))
    print('Nucleotide genes written:', found)
    print('Missing roots:', len(missing))
    if missing:
        print('First missing roots:', ', '.join(sorted(missing)[:20]))
        raise SystemExit('ERROR: Some final gene roots were not found in the nucleotide gene FASTA.')
    if args.expected is not None and found != args.expected:
        raise SystemExit(f'ERROR: Expected {args.expected} genes, wrote {found}.')


if __name__ == '__main__':
    main()
