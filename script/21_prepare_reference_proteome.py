#!/usr/bin/env python3
import argparse, csv
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    d={}; h=None; seq=[]
    with open(path) as fh:
        for line in fh:
            line=line.rstrip('\n')
            if not line: continue
            if line.startswith('>'):
                if h is not None: d[h.split()[0]]=(''.join(seq), h)
                h=line[1:]; seq=[]
            else: seq.append(line.strip())
    if h is not None: d[h.split()[0]]=(''.join(seq), h)
    return d


def attrs(s):
    d={}
    for x in s.split(';'):
        if '=' in x:
            k,v=x.split('=',1); d[k]=v
    return d


def wrap(s,n=80): return '\n'.join(s[i:i+n] for i in range(0,len(s),n))


def main():
    p=argparse.ArgumentParser(description='Prepare one representative protein per NCBI locus_tag while optionally excluding organelle contigs.')
    p.add_argument('--proteins', required=True)
    p.add_argument('--gff', required=True)
    p.add_argument('--out-fasta', required=True)
    p.add_argument('--out-map', required=True)
    p.add_argument('--exclude-contigs', default='', help='Comma-separated contig accessions to exclude.')
    p.add_argument('--id-mode', choices=['locus','protein'], default='locus', help='Use locus_tag or source protein_id as output FASTA ID.')
    args=p.parse_args()
    exclude={x for x in args.exclude_contigs.split(',') if x}
    seqs=read_fasta(args.proteins)

    locus_to_pids=defaultdict(set); pid_to_contig={}
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith('#'): continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)!=9 or parts[2]!='CDS' or parts[0] in exclude: continue
            a=attrs(parts[8]); locus=a.get('locus_tag'); pid=a.get('protein_id')
            if locus and pid:
                locus_to_pids[locus].add(pid); pid_to_contig[pid]=parts[0]

    reps=[]; missing=[]
    for locus,pids in locus_to_pids.items():
        available=[pid for pid in pids if pid in seqs]
        if not available:
            missing.append(locus); continue
        pid=max(available, key=lambda x: len(seqs[x][0]))
        seq=seqs[pid][0]
        reps.append((locus,pid,seq,pid_to_contig.get(pid,''),len(pids)))
    reps.sort(key=lambda x:x[0])

    Path(args.out_fasta).parent.mkdir(parents=True,exist_ok=True)
    with open(args.out_fasta,'w') as out:
        for locus,pid,seq,contig,nmodels in reps:
            out_id=locus if args.id_mode=='locus' else pid
            out.write(f'>{out_id} locus_tag={locus} representative_protein_id={pid}\n{wrap(seq)}\n')
    with open(args.out_map,'w',newline='') as out:
        w=csv.writer(out,delimiter='\t')
        w.writerow(['locus_tag','representative_protein_id','source_contig','protein_length_aa','protein_models_for_locus'])
        for locus,pid,seq,contig,nmodels in reps:
            w.writerow([locus,pid,contig,len(seq),nmodels])

    print('Input protein FASTA records:', len(seqs))
    print('Excluded contigs:', ','.join(sorted(exclude)) if exclude else 'none')
    print('Protein-coding loci retained:', len(locus_to_pids))
    print('Loci with >1 protein model:', sum(len(v)>1 for v in locus_to_pids.values()))
    print('Representative proteins written:', len(reps))
    print('Loci missing FASTA sequence:', len(missing))

if __name__=='__main__':
    main()
