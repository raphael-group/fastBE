import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

def draw_graph(adjacency_list_file, ax=None):
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
    nx.draw(T, pos, with_labels=True, arrows=True, ax=ax)

def main():
    parser = argparse.ArgumentParser(description="Draw a graph from an adjacency list")
    parser.add_argument("adjacency_lists", help="The file containing the adjacency list", nargs='+')

    args = parser.parse_args()

    fig, axes = plt.subplots(nrows=1, ncols=len(args.adjacency_lists))
    for i, adjacency_list_file in enumerate(args.adjacency_lists):
        axes[i].set_title(adjacency_list_file)
        draw_graph(adjacency_list_file, axes[i])

    plt.show()
if __name__ == "__main__":
    main()

