import pandas as pd
import networkx as nx

edge_file = "hic_intercontig_contacts_all_primary_pairs.tsv"

gexf_out = "hic_contig_network_all_primary_pairs.gexf"
graphml_out = "hic_contig_network_all_primary_pairs.graphml"

edges = pd.read_csv(edge_file, sep="\t", dtype=str)

expected = {"contig_A", "contig_B", "HiC_contact_pairs"}
if not expected.issubset(set(edges.columns)):
    edges = pd.read_csv(
        edge_file,
        sep="\t",
        header=None,
        names=["contig_A", "contig_B", "HiC_contact_pairs"],
        dtype=str
    )

edges = edges[edges["contig_A"] != "contig_A"].copy()
edges["HiC_contact_pairs"] = pd.to_numeric(
    edges["HiC_contact_pairs"],
    errors="coerce"
).fillna(0).astype(int)

G = nx.Graph()

for _, row in edges.iterrows():
    a = str(row["contig_A"])
    b = str(row["contig_B"])
    weight = int(row["HiC_contact_pairs"])

    if a == b:
        continue

    G.add_edge(
        a,
        b,
        weight=weight,
        HiC_contact_pairs=weight
    )

for node in G.nodes():
    edge_weights = [G[node][nbr]["weight"] for nbr in G.neighbors(node)]

    G.nodes[node]["connected_partners"] = len(edge_weights)
    G.nodes[node]["total_intercontig_HiC_pairs"] = sum(edge_weights)
    G.nodes[node]["strongest_single_contact_pairs"] = max(edge_weights) if edge_weights else 0

print(f"Nodes: {G.number_of_nodes()}")
print(f"Edges: {G.number_of_edges()}")

nx.write_gexf(G, gexf_out)
nx.write_graphml(G, graphml_out)

print(f"Wrote: {gexf_out}")
print(f"Wrote: {graphml_out}")
