import pandas as pd
import networkx as nx
import json 
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("adjacency_list", type=str, help="Adjacency list defining the tree.")
    parser.add_argument("params", type=str, help="Pairtree parameters as a JSON file.")
    parser.add_argument("-o", "--output", type=str, required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    with open(args.params, "r") as f:
        params = json.load(f)

    T = nx.read_adjlist(args.adjacency_list, nodetype=int, create_using=nx.DiGraph)
    parents = [next(T.predecessors(i)) for i in range(1, len(T.nodes))]
    params['structures'] = [parents]

    with open(args.output, "w") as f:
        json.dump(params, f)
