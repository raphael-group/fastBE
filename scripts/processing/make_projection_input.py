import numpy as np
import pandas as pd
import networkx as nx
import argparse

# input files: examples/sim_tree.txt examples/sim_obs_frequency_matrix.txt
# output files: examples/sim_obs_projection.txt

def parse_args():
    parser = argparse.ArgumentParser(description='Make input for PPM projection algorithm.')
    parser.add_argument('tree', type=str, help='Tree file')
    parser.add_argument('frequency_matrix', type=str, help='Frequency matrix file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
    frequency_matrix = np.loadtxt(args.frequency_matrix)

    if len(frequency_matrix.shape) == 1:
        frequency_matrix = frequency_matrix.reshape(1, -1)

    m, n = frequency_matrix.shape

    print(f"{n} {m}")
    for j in range(n):
        for i in range(m):
            print(frequency_matrix[i,j], end=" ")
    print()

    for i in range(n):
        print("1.0", end=" ")
    print()

    root = [node for node in tree.nodes if tree.in_degree(node) == 0][0]
    print(root)
    for i in range(n):
        deg = tree.degree(i)
        print(deg, end=" ")

    print()
    for i in range(n):
        neighbors = list(tree.successors(i)) + list(tree.predecessors(i))
        print(" ".join(map(str, neighbors)))
    print(0)

