import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout

def draw_graph(adjacency_list_file):
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
    nx.draw(T, pos, with_labels=True, arrows=True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Draw a graph from an adjacency list")
    parser.add_argument("adjacency_list_file", help="The file containing the adjacency list")

    args = parser.parse_args()

    draw_graph(args.adjacency_list_file)

if __name__ == "__main__":
    main()

