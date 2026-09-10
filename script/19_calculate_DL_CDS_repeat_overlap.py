#!/usr/bin/env python3
import argparse, csv, heapq, os, sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


def parse_attributes(s):
    d = {}
    for item in s.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            d[k] = v
    return d


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def coarse_class(c):
    if c == 'Unknown': return 'Unknown'
    if c.startswith('LTR/Copia'): return 'LTR_Copia'
    if c.startswith('LTR/Gypsy'): return 'LTR_Gypsy'
    if c.startswith('LTR/Ngaro'): return 'LTR_Ngaro'
    if c.startswith('LTR/'): return 'LTR_other'
    if c.startswith('LINE/'): return 'LINE'
    if c.startswith('DNA/PIF-Harb'): return 'DNA_PIF_Harbinger'
    if c.startswith('DNA/'): return 'DNA_other'
    if c == 'Simple_repeat': return 'Simple_repeat'
    if c == 'Low_complexity': return 'Low_complexity'
    if c.startswith('Satellite'): return 'Satellite'
    return 'Other_repeat'


INTERSPERSED = {'Unknown','LTR_Copia','LTR_Gypsy','LTR_Ngaro','LTR_other','LINE','DNA_PIF_Harbinger','DNA_other','Other_repeat'}
BIN_ORDER = ['0','>0-10','>10-25','>25-50','>50-75','>75-100']


def overlap_bin(v):
    if v == 0: return '0'
    if v <= 10: return '>0-10'
    if v <= 25: return '>10-25'
    if v <= 50: return '>25-50'
    if v <= 75: return '>50-75'
    return '>75-100'


def resolve_repeat_hits(hits):
    events = defaultdict(lambda: {'start': [], 'end': []})
    for i, h in enumerate(hits):
        s, e, score, name, cls = h
        events[s]['start'].append(i)
        events[e]['end'].append(i)
    active, heap, segments = set(), [], []
    prev = None
    for pos in sorted(events):
        if prev is not None and pos > prev and active:
            while heap and heap[0][1] not in active:
                heapq.heappop(heap)
            if heap:
                i = heap[0][1]
                s, e, score, name, cls = hits[i]
                seg = (prev, pos, cls, name, score)
                if segments and segments[-1][1] == seg[0] and segments[-1][2:4] == seg[2:4]:
                    old = segments[-1]
                    segments[-1] = (old[0], seg[1], old[2], old[3], max(old[4], seg[4]))
                else:
                    segments.append(seg)
        for i in events[pos]['end']:
            active.discard(i)
        for i in events[pos]['start']:
            active.add(i)
            heapq.heappush(heap, (-hits[i][2], i))
        prev = pos
    return segments


def main():
    p = argparse.ArgumentParser(description='Calculate RepeatMasker overlap across CDS bases of final Deer Lake representative transcripts.')
    p.add_argument('--map', required=True, help='Representative proteome map TSV.')
    p.add_argument('--gff', required=True, help='BRAKER GFF3.')
    p.add_argument('--repeatmasker-out', required=True, help='RepeatMasker .out file.')
    p.add_argument('--outdir', required=True)
    p.add_argument('--candidate-threshold', type=float, default=50.0)
    args = p.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    with open(args.map) as fh:
        map_rows = list(csv.DictReader(fh, delimiter='\t'))
    selected = {r['gene_id'] for r in map_rows}

    cds = defaultdict(list); contigs = defaultdict(set)
    with open(args.gff) as fh:
        for line in fh:
            if line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9 or parts[2] != 'CDS': continue
            a = parse_attributes(parts[8])
            tx = a.get('transcript_id') or a.get('Parent')
            if tx not in selected: continue
            cds[tx].append((int(parts[3])-1, int(parts[4])))
            contigs[tx].add(parts[0])
    missing = selected - set(cds)
    if missing:
        sys.exit(f'ERROR: {len(missing)} representative transcripts have no CDS in GFF.')
    if any(len(v) != 1 for v in contigs.values()):
        sys.exit('ERROR: at least one transcript occurs on multiple contigs.')
    for tx in cds:
        cds[tx] = merge_intervals(cds[tx])

    hits_by_contig = defaultdict(list)
    with open(args.repeatmasker_out) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 15: continue
            try:
                score = int(parts[0]); start = int(parts[5])-1; end = int(parts[6])
            except ValueError:
                continue
            hits_by_contig[parts[4]].append((start, end, score, parts[9], parts[10]))

    resolved = {}; starts = {}
    for c, hits in hits_by_contig.items():
        x = resolve_repeat_hits(hits)
        resolved[c] = x
        starts[c] = [s[0] for s in x]

    results, breakdown = [], []
    for row in map_rows:
        root, tx, contig = row['gene_root'], row['gene_id'], row['contig_id']
        gff_contig = next(iter(contigs[tx]))
        if contig != gff_contig:
            sys.exit(f'ERROR: contig mismatch for {tx}: map={contig}, GFF={gff_contig}')
        cds_len = sum(e-s for s,e in cds[tx])
        cat_bp, raw_bp = defaultdict(int), defaultdict(int)
        segs, sstarts = resolved.get(contig, []), starts.get(contig, [])
        for cs, ce in cds[tx]:
            if not segs: continue
            i = max(0, bisect_right(sstarts, cs)-1)
            while i < len(segs) and segs[i][1] <= cs: i += 1
            while i < len(segs) and segs[i][0] < ce:
                rs, re, cls, name, score = segs[i]
                ov = min(ce,re) - max(cs,rs)
                if ov > 0:
                    cat_bp[coarse_class(cls)] += ov
                    raw_bp[cls] += ov
                i += 1
        total_bp = sum(cat_bp.values())
        inter_bp = sum(v for k,v in cat_bp.items() if k in INTERSPERSED)
        total_pct = 100*total_bp/cds_len if cds_len else 0
        inter_pct = 100*inter_bp/cds_len if cds_len else 0
        dom = max(raw_bp, key=raw_bp.get) if raw_bp else ''
        out = {
            'gene_root':root,'transcript_id':tx,'contig_id':contig,'protein_length_aa':row['protein_length_aa'],
            'diatom_Average_TPM':row['diatom_Average_TPM'],'functional_annotation':row['functional_annotation'],
            'CDS_length_bp':cds_len,'repeat_overlap_bp':total_bp,'repeat_overlap_percent':round(total_pct,4),
            'interspersed_repeat_overlap_bp':inter_bp,'interspersed_repeat_overlap_percent':round(inter_pct,4),
            'LTR_Copia_bp':cat_bp['LTR_Copia'],'LTR_Gypsy_bp':cat_bp['LTR_Gypsy'],'LTR_Ngaro_bp':cat_bp['LTR_Ngaro'],
            'LTR_other_bp':cat_bp['LTR_other'],'LINE_bp':cat_bp['LINE'],'DNA_PIF_Harbinger_bp':cat_bp['DNA_PIF_Harbinger'],
            'DNA_other_bp':cat_bp['DNA_other'],'Unknown_bp':cat_bp['Unknown'],'Simple_repeat_bp':cat_bp['Simple_repeat'],
            'Low_complexity_bp':cat_bp['Low_complexity'],'Satellite_bp':cat_bp['Satellite'],'Other_repeat_bp':cat_bp['Other_repeat'],
            'dominant_repeat_class':dom,'dominant_repeat_class_bp':raw_bp.get(dom,0)
        }
        results.append(out)
        for cls,bp in sorted(raw_bp.items()):
            breakdown.append({'gene_root':root,'transcript_id':tx,'repeat_class':cls,'overlap_bp':bp})

    fields = list(results[0].keys())
    main_out = outdir/'DL_representative_proteins_CDS_repeat_overlap.tsv'
    with open(main_out,'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(results)
    with open(outdir/'DL_representative_CDS_repeat_class_breakdown.tsv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['gene_root','transcript_id','repeat_class','overlap_bp'],delimiter='\t'); w.writeheader(); w.writerows(breakdown)

    summary=[]
    for metric in ['repeat_overlap_percent','interspersed_repeat_overlap_percent']:
        counts=defaultdict(int)
        for r in results: counts[overlap_bin(float(r[metric]))]+=1
        for b in BIN_ORDER:
            summary.append({'metric':metric,'bin':b,'gene_count':counts[b],'percentage_of_total':round(100*counts[b]/len(results),4)})
    with open(outdir/'DL_CDS_repeat_overlap_bins.tsv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['metric','bin','gene_count','percentage_of_total'],delimiter='\t'); w.writeheader(); w.writerows(summary)

    high=[r for r in results if float(r['interspersed_repeat_overlap_percent']) >= args.candidate_threshold]
    high.sort(key=lambda r: float(r['interspersed_repeat_overlap_percent']), reverse=True)
    with open(outdir/'DL_CDS_repeat_overlap_candidates_ge50.tsv','w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(high)

    print('Representative proteins:', len(results))
    print('Genes with any interspersed-repeat CDS overlap:', sum(float(r['interspersed_repeat_overlap_percent'])>0 for r in results))
    print(f'Genes with >={args.candidate_threshold:g}% interspersed-repeat CDS overlap:', len(high))
    print('Main output:', main_out)

if __name__ == '__main__':
    main()
