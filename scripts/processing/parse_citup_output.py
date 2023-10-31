import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Parse CITUP output')
    parser.add_argument('results', type=str, help='HDF file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    tree_ids = pd.read_hdf(args.results, '/results/tree_id')
    lh = -1 * pd.read_hdf(args.results, '/results/likelihood')

    best_tree = tree_ids[lh.idxmax()]
    edge_list = pd.read_hdf(args.results, "/trees/" + str(best_tree) + "/adjacency_list").to_numpy()
    adjacency_list = [[] for _ in range(len(edge_list) + 1)]
    for edge in edge_list:
        adjacency_list[edge[0]].append(edge[1])

    for i in range(len(adjacency_list)):
        print(str(i) + " " + " ".join(map(str, adjacency_list[i])))
