#!/usr/bin/env python3
import argparse, csv, re, sys
from pathlib import Path


def read_fasta(path):
    records = {}
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records[header.split()[0]] = ''.join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if header is not None:
        records[header.split()[0]] = ''.join(seq)
    return records


def wrap(seq, width=80):
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))


def gene_root(transcript_id):
    return re.sub(r'\.t\d+$', '', transcript_id)


def has_tpm(value):
    return value is not None and value.strip() != ''


def low_information_annotation(value):
    x = (value or '').strip().lower()
    return x in {'', 'unannotated', 'uncharacterized protein'}


def main():
    p = argparse.ArgumentParser(description='Create one representative Deer Lake nuclear protein per BRAKER gene.')
    p.add_argument('--gene-table', required=True, help='TSV with gene_id, contig_id, diatom_compartment, protein/annotation/TPM fields.')
    p.add_argument('--proteins', required=True, help='BRAKER protein FASTA.')
    p.add_argument('--out-fasta', required=True)
    p.add_argument('--out-map', required=True)
    p.add_argument('--antifam-gene-roots', default='g10893,g11404', help='Comma-separated gene roots to exclude.')
    p.add_argument('--min-aa', type=int, default=50, help='Remove proteins shorter than this only when unannotated and without TPM.')
    args = p.parse_args()

    proteins = read_fasta(args.proteins)
    antifam = {x for x in args.antifam_gene_roots.split(',') if x}

    with open(args.gene_table) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))

    required = {'gene_id', 'contig_id', 'diatom_compartment', 'functional_annotation', 'diatom_Average_TPM'}
    missing = required - set(rows[0].keys()) if rows else required
    if missing:
        sys.exit(f'ERROR: missing columns: {sorted(missing)}')

    nuclear = [r for r in rows if r['diatom_compartment'].strip().lower() == 'nuclear']
    candidates = {}
    for r in nuclear:
        tx = r['gene_id'].strip()
        root = gene_root(tx)
        if root in antifam:
            continue
        seq = proteins.get(tx)
        if not seq:
            continue
        if '*' in seq:
            continue
        prev = candidates.get(root)
        if prev is None or len(seq) > len(prev['seq']):
            candidates[root] = {'row': r, 'tx': tx, 'seq': seq}

    kept = []
    removed_short = []
    for root, d in candidates.items():
        r, seq = d['row'], d['seq']
        if len(seq) < args.min_aa and not has_tpm(r['diatom_Average_TPM']) and low_information_annotation(r['functional_annotation']):
            removed_short.append((root, d['tx'], len(seq)))
            continue
        kept.append((root, d))

    kept.sort(key=lambda x: x[0])
    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_map).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_fasta, 'w') as out:
        for root, d in kept:
            out.write(f'>{root} representative_isoform={d["tx"]}\n{wrap(d["seq"])}\n')

    fields = ['gene_root', 'gene_id', 'contig_id', 'diatom_compartment', 'protein_length_aa', 'diatom_Average_TPM', 'functional_annotation']
    with open(args.out_map, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for root, d in kept:
            r = d['row']
            w.writerow({
                'gene_root': root,
                'gene_id': d['tx'],
                'contig_id': r['contig_id'],
                'diatom_compartment': r['diatom_compartment'],
                'protein_length_aa': len(d['seq']),
                'diatom_Average_TPM': r['diatom_Average_TPM'],
                'functional_annotation': r['functional_annotation'],
            })

    print('Input table rows:', len(rows))
    print('Nuclear isoform rows:', len(nuclear))
    print('AntiFam gene roots excluded:', len(antifam))
    print('Representative genes before short-protein filter:', len(candidates))
    print('Short unsupported proteins removed:', len(removed_short))
    print('Final representative proteins:', len(kept))
    if removed_short:
        print('Removed short proteins:')
        for x in removed_short:
            print('\t'.join(map(str, x)))


if __name__ == '__main__':
    main()
