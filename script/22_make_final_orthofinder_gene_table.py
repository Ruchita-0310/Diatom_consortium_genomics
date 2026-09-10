#!/usr/bin/env python3
import argparse, csv, sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description='Create the final Deer Lake gene table with TPM and four OrthoFinder orthogroup-sharing yes/no fields.')
    p.add_argument('--orthogroups', required=True, help='OrthoFinder Orthogroups.tsv')
    p.add_argument('--dl-map', required=True, help='Final Deer Lake representative proteome map TSV')
    p.add_argument('--metadata', required=True, help='Existing Deer Lake metadata/TPM table keyed by gene_id')
    p.add_argument('--out', required=True)
    p.add_argument('--expected-genes', type=int, default=None)
    args=p.parse_args()

    with open(args.metadata) as fh:
        metadata={r['gene_id']:r for r in csv.DictReader(fh,delimiter='\t')}
    with open(args.dl_map) as fh:
        final_genes=[{'gene_root':r['gene_root'],'gene_id':r['gene_id']} for r in csv.DictReader(fh,delimiter='\t')]
    if args.expected_genes is not None and len(final_genes)!=args.expected_genes:
        sys.exit(f'ERROR: expected {args.expected_genes} Deer Lake genes, found {len(final_genes)}')

    presence={g['gene_root']:{'TP':'no','PT':'no','NI':'no','SR':'no'} for g in final_genes}
    with open(args.orthogroups) as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        required={'DeerLake_Nitzschia','Nitzschia_inconspicua','Phaeodactylum_tricornutum','Seminavis_robusta','Thalassiosira_pseudonana'}
        if not required.issubset(reader.fieldnames or []):
            sys.exit('ERROR: OrthoFinder species columns do not match the expected five-species run.')
        for row in reader:
            dl=row['DeerLake_Nitzschia'].strip()
            if not dl: continue
            state={
                'TP':'yes' if row['Thalassiosira_pseudonana'].strip() else 'no',
                'PT':'yes' if row['Phaeodactylum_tricornutum'].strip() else 'no',
                'NI':'yes' if row['Nitzschia_inconspicua'].strip() else 'no',
                'SR':'yes' if row['Seminavis_robusta'].strip() else 'no',
            }
            for root in [x.strip() for x in dl.split(',') if x.strip()]:
                if root in presence: presence[root]=state.copy()

    fields=['gene_id','contig_id','diatom_compartment','diatom_gene_length_bp','functional_annotation','diatom_Average_TPM',
            'present_in_Thalassiosira_pseudonana','present_in_Phaeodactylum_tricornutum','present_in_Nitzschia_inconspicua','present_in_Seminavis_robusta']
    rows=[]; missing=[]
    for g in final_genes:
        root,tx=g['gene_root'],g['gene_id']
        if tx not in metadata:
            missing.append(tx); continue
        m=metadata[tx]; q=presence[root]
        rows.append({
            'gene_id':tx,'contig_id':m['contig_id'],'diatom_compartment':m['diatom_compartment'],
            'diatom_gene_length_bp':m['diatom_gene_length_bp'],'functional_annotation':m['functional_annotation'],
            'diatom_Average_TPM':m['diatom_Average_TPM'],
            'present_in_Thalassiosira_pseudonana':q['TP'],'present_in_Phaeodactylum_tricornutum':q['PT'],
            'present_in_Nitzschia_inconspicua':q['NI'],'present_in_Seminavis_robusta':q['SR']})
    if missing:
        sys.exit(f'ERROR: missing metadata for {len(missing)} representative transcripts; examples: {missing[:10]}')

    Path(args.out).parent.mkdir(parents=True,exist_ok=True)
    with open(args.out,'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(rows)

    print('Final Deer Lake genes:', len(final_genes))
    print('Final table rows:', len(rows))
    for label,col in [('TP','present_in_Thalassiosira_pseudonana'),('PT','present_in_Phaeodactylum_tricornutum'),('NI','present_in_Nitzschia_inconspicua'),('SR','present_in_Seminavis_robusta')]:
        print(f'{label} yes:', sum(r[col]=='yes' for r in rows))
    print('Saved:', args.out)

if __name__=='__main__':
    main()
