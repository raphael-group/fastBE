import random
import pandas as pd
import numpy as np
import json
import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from matplotlib.patches import Rectangle

def draw_graph(adjacency_list_file, usage_matrix, ssm_df, clusters, ax=None):
    T = nx.DiGraph()

    with open(adjacency_list_file, 'r') as f:
        for line in f:
            # Ignore lines starting with '#'
            if line[0] == '#':
                continue
            
            nodes = list(map(int, line.split()))
            parent_node = nodes[0]
            
            for node in nodes[1:]:
                T.add_edge(parent_node, node)

    random.seed(0)
    color_map = {}
    for node in sorted(T.nodes()):
        color = random.choice([2,3,4])
        color_map[node] = color

    print("digraph G {")
    print("fontname = \"helvetica\";")
    print("edge[arrowhead=\"vee\"];")
    print("node[colorscheme=gnbu9]")
    for node in T.nodes():
        if node == 0:
            print("{} [shape=box style=\"rounded,filled\" fillcolor=5 penwidth=0 label=<<B><I>Normal Clone</I></B>>]".format(node))
            continue

        mutations = clusters[node - 1]
        mutation_records = ssm_df.loc[mutations]
        record = mutation_records.name.tolist()

        html_record_str = f"<I><B>Clone {node}:</B></I><BR/>"
        for i in range(len(record)):
            html_record_str += "{}".format(record[i].split("_")[0])
            if i != len(record) - 1:
                html_record_str += "<BR/>"
        color = color_map[node]
        print("{} [shape=box style=\"rounded,filled\" fillcolor={} penwidth=0 label=<{}>]".format(node, color, html_record_str))

    for edge in T.edges():
        print("{} -> {}".format(edge[0], edge[1]))
    print("}")

def main():
    parser = argparse.ArgumentParser(description="Draw a graph from an adjacency list")
    parser.add_argument("adjacency_list", help="The file containing the adjacency list")
    parser.add_argument("usage_matrix", help="Usage matrix")
    parser.add_argument("ssm_file", help="Simple somatic mutation file")
    parser.add_argument("params", help="Parameters file (JSON)")
    args = parser.parse_args()

    usage_matrix = np.loadtxt(args.usage_matrix, dtype=float, delimiter=' ')
    ssm_df = pd.read_csv(args.ssm_file, sep='\t', header=0, index_col=0)
    with open(args.params, 'r') as f:
        params = json.load(f)
    draw_graph(args.adjacency_list, usage_matrix, ssm_df, params['clusters'])

if __name__ == "__main__":
    main()

