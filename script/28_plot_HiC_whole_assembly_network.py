#!/usr/bin/env python3
"""Plot all polished-assembly contigs and the retained high-confidence Hi-C network."""

import argparse
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DIATOM = '#0072B2'
BACTERIAL = '#E69F00'
PLASTID = '#009E73'
MITO = '#CC79A7'
UNKNOWN = '#999999'
DIATOM_EDGE = '#56B4E9'
BACTERIAL_EDGE = '#E69F00'
MIXED_EDGE = '#4D4D4D'
PLASTID_CONTIGS = {'contig_1443','contig_4315'}
MITO_CONTIGS = {'contig_5628','contig_1647'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--edges', required=True)
    p.add_argument('--types', required=True)
    p.add_argument('--out-prefix', default='HiC_whole_assembly_network')
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    edges = pd.read_csv(args.edges, sep='\t')
    types = pd.read_csv(args.types, sep='\t')
    types['contig_id'] = types['contig_id'].astype(str).str.strip()
    types['contig_type'] = types['contig_type'].astype(str).str.strip().str.lower()

    G = nx.Graph()
    for _, r in types.iterrows():
        G.add_node(r['contig_id'], contig_type=r['contig_type'])
    for _, r in edges.iterrows():
        G.add_edge(str(r['node1']), str(r['node2']), weight=float(r['read_pair_count']))

    connected = [n for n in G if G.degree(n) > 0]
    isolated = [n for n in G if G.degree(n) == 0]
    diatom = [n for n in G if G.nodes[n].get('contig_type') == 'diatom']
    bacterial = [n for n in G if G.nodes[n].get('contig_type') == 'bacterial']
    unknown = [n for n in G if G.nodes[n].get('contig_type') not in {'diatom','bacterial'}]
    plastid = [n for n in PLASTID_CONTIGS if n in G]
    mito = [n for n in MITO_CONTIGS if n in G]

    print('Total assembly contigs:', G.number_of_nodes())
    print('Hi-C edges:', G.number_of_edges())
    print('Connected contigs:', len(connected))
    print('Isolated contigs:', len(isolated))
    print('Diatom contigs:', len(diatom))
    print('Bacterial contigs:', len(bacterial))
    print('Unknown contigs:', len(unknown))
    print('Plastid contigs found:', plastid)
    print('Candidate mitochondrial contigs found:', mito)
    for n in sorted(set(plastid + mito)):
        print(n, '| type:', G.nodes[n].get('contig_type'), '| degree:', G.degree(n), '| weighted degree:', G.degree(n, weight='weight'))

    # Layout only connected contigs.
    CG = G.subgraph(connected).copy()
    pos = nx.spring_layout(CG, seed=args.seed, weight='weight', k=0.25, iterations=500)
    pos = nx.rescale_layout_dict(pos, scale=1.0)

    # Place isolates in concentric rings to show all assembly contigs without implying contacts.
    if isolated:
        rng = np.random.default_rng(args.seed)
        rings = np.array_split(np.array(isolated, dtype=object), 5)
        radii = np.linspace(1.6, 2.8, 5)
        for ring_nodes, radius in zip(rings, radii):
            ring_nodes = list(ring_nodes)
            rng.shuffle(ring_nodes)
            angles = np.linspace(0, 2*np.pi, len(ring_nodes), endpoint=False) + rng.uniform(0, 2*np.pi)
            for n, angle in zip(ring_nodes, angles):
                pos[n] = (radius*np.cos(angle), radius*np.sin(angle))

    conn_d = [n for n in connected if G.nodes[n].get('contig_type') == 'diatom']
    conn_b = [n for n in connected if G.nodes[n].get('contig_type') == 'bacterial']
    conn_u = [n for n in connected if G.nodes[n].get('contig_type') not in {'diatom','bacterial'}]
    iso_d = [n for n in isolated if G.nodes[n].get('contig_type') == 'diatom']
    iso_b = [n for n in isolated if G.nodes[n].get('contig_type') == 'bacterial']
    iso_u = [n for n in isolated if G.nodes[n].get('contig_type') not in {'diatom','bacterial'}]

    dd, bb, mixed = [], [], []
    for u,v,d in G.edges(data=True):
        tu, tv = G.nodes[u].get('contig_type','unknown'), G.nodes[v].get('contig_type','unknown')
        item = (u,v,float(d.get('weight',1)))
        if tu == tv == 'diatom': dd.append(item)
        elif tu == tv == 'bacterial': bb.append(item)
        else: mixed.append(item)
    maxw = max([w for _,_,w in dd+bb+mixed], default=1)
    def width(w): return 0.25 + 2.0*np.log10(w+1)/np.log10(maxw+1)

    fig, ax = plt.subplots(figsize=(12,12))
    # Isolated nodes first.
    for nodes, color in [(iso_d,DIATOM),(iso_b,BACTERIAL),(iso_u,UNKNOWN)]:
        if nodes: nx.draw_networkx_nodes(G,pos,nodelist=nodes,node_size=7,node_color=color,alpha=0.30,linewidths=0,ax=ax)
    # Edges.
    for items, color, alpha in [(dd,DIATOM_EDGE,0.32),(bb,BACTERIAL_EDGE,0.30),(mixed,MIXED_EDGE,0.65)]:
        if items:
            nx.draw_networkx_edges(G,pos,edgelist=[(u,v) for u,v,_ in items],width=[width(w) for _,_,w in items],edge_color=color,alpha=alpha,ax=ax)
    # Connected nodes.
    for nodes, color in [(conn_d,DIATOM),(conn_b,BACTERIAL),(conn_u,UNKNOWN)]:
        if nodes: nx.draw_networkx_nodes(G,pos,nodelist=nodes,node_size=22,node_color=color,alpha=0.90,linewidths=0.2,edgecolors='black',ax=ax)
    # Organelle stars.
    if plastid: nx.draw_networkx_nodes(G,pos,nodelist=plastid,node_size=220,node_color=PLASTID,node_shape='*',edgecolors='black',linewidths=1.0,alpha=1.0,ax=ax)
    if mito: nx.draw_networkx_nodes(G,pos,nodelist=mito,node_size=220,node_color=MITO,node_shape='*',edgecolors='black',linewidths=1.0,alpha=1.0,ax=ax)
    labels = {n:n for n in plastid+mito}
    if labels: nx.draw_networkx_labels(G,pos,labels=labels,font_size=8,font_weight='bold',verticalalignment='bottom',horizontalalignment='left',ax=ax)

    legend = [
        Line2D([0],[0],marker='o',linestyle='None',markerfacecolor=DIATOM,markeredgecolor='black',markeredgewidth=0.3,markersize=8,label=f'Diatom contigs (n = {len(diatom):,})'),
        Line2D([0],[0],marker='o',linestyle='None',markerfacecolor=BACTERIAL,markeredgecolor='black',markeredgewidth=0.3,markersize=8,label=f'Bacterial contigs (n = {len(bacterial):,})'),
        Line2D([0],[0],marker='*',linestyle='None',markerfacecolor=PLASTID,markeredgecolor='black',markersize=13,label='Plastid genome associated contigs'),
        Line2D([0],[0],marker='*',linestyle='None',markerfacecolor=MITO,markeredgecolor='black',markersize=13,label='Candidate mitochondrial contigs'),
        Line2D([0],[0],color=MIXED_EDGE,lw=2,label='Diatom-bacterial Hi-C contact'),
    ]
    ax.legend(handles=legend,loc='lower left',frameon=False,fontsize=10)
    ax.text(0.02,0.98,f'Assembly contigs = {G.number_of_nodes():,}\nHi-C connected contigs = {len(connected):,}\nHi-C edges = {G.number_of_edges():,}\nMinimum edge support = 2 read pairs',transform=ax.transAxes,ha='left',va='top',fontsize=10)
    ax.text(-0.01,1.01,'B',transform=ax.transAxes,fontsize=22,fontweight='bold',ha='left',va='bottom')
    ax.set_axis_off(); ax.set_aspect('equal'); plt.tight_layout()
    for ext in ['png','pdf','svg']:
        kwargs={'bbox_inches':'tight'}
        if ext=='png': kwargs['dpi']=600
        plt.savefig(f'{args.out_prefix}.{ext}', **kwargs)
    plt.show()


if __name__ == '__main__':
    main()
