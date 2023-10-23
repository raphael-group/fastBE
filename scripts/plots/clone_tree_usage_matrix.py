import numpy as np
import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
from matplotlib.patches import Rectangle

def draw_graph(adjacency_list_file, usage_matrix, ax=None):
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
    
    pos = graphviz_layout(T, prog='dot')

    # Enhanced aesthetics for nodes and edges
    node_size = 800
    node_color = 'skyblue'
    edge_color = 'gray'
    arrow_size = 20
    edge_alpha = 0.6

    labels = {}
    for node in T.nodes():
        labels[node] = node 

    nx.draw(T, pos, with_labels=False, arrows=True, node_size=node_size, node_color=node_color,
            edge_color=edge_color, alpha=edge_alpha, ax=ax, arrowsize=arrow_size)
    nx.draw_networkx_labels(T, pos, labels, font_size=12, ax=ax)

    entry_width = 10
    box_height = 10
    entry_spacing = 3

    # Adjusted vector box aesthetics
    box_edge_color = 'black'
    box_face_color = 'lightgray'
    text_color = 'black'

    for node in T.nodes():
        if T.out_degree(node) == 0:
            x_offset, y_offset = -25, -20  # Adjusted offset
        else:
            x_offset, y_offset = -60, 0  # Adjusted offset

        vector_length = len(usage_matrix[:, node])
        
        box_width = vector_length * entry_width + (vector_length - 1) * entry_spacing
        
        rect_x = pos[node][0] + x_offset
        rect_y = pos[node][1] - box_height / 2 + y_offset
        rect = Rectangle((rect_x, rect_y), box_width, box_height, facecolor=box_face_color, edgecolor=box_edge_color, linewidth=1.2)
        ax.add_patch(rect)
        
        for idx, entry in enumerate(usage_matrix[:, node]):
            x = rect_x + idx * (entry_width + entry_spacing) + entry_width / 2
            y = rect_y + box_height / 2
            ax.text(x, y, f"{entry:.2f}", ha='center', va='center', fontsize=10, color=text_color)

def main():
    parser = argparse.ArgumentParser(description="Draw a graph from an adjacency list")
    parser.add_argument("adjacency_list", help="The file containing the adjacency list")
    parser.add_argument("usage_matrix", help="Usage matrix")
    args = parser.parse_args()

    usage_matrix = np.loadtxt(args.usage_matrix, dtype=float, delimiter=' ')

    # entrywise log and dot with itself
    # replace Nan with 0
    shannon_diversity = -np.log(usage_matrix) * usage_matrix 
    shannon_diversity[np.isnan(shannon_diversity)] = 0
    shannon_diversity = np.sum(shannon_diversity, axis=1)
    print(shannon_diversity)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    draw_graph(args.adjacency_list, usage_matrix, ax)
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

