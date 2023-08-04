import numpy as np

import argparse

def convert_parents_to_adjmatrix(parents):
  K = len(parents) + 1
  adjm = np.eye(K)
  adjm[parents,np.arange(1, K)] = 1
  return adjm

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse pairtree output')
    parser.add_argument('result', type=str, help='Result file in NPZ format.')
    parser.add_argument('--output', type=str, help='Output adjacency list.')
    args = parser.parse_args()

    result = np.load(args.result)
    best_tree_adj = convert_parents_to_adjmatrix(result['struct'][0])

    adjacency_list = [[] for _ in range(best_tree_adj.shape[0])]
    for i in range(best_tree_adj.shape[0]):
        for j in range(best_tree_adj.shape[1]):
            if best_tree_adj[i, j] == 1:
                adjacency_list[i].append(j)

    # write adjacency list to NX format
    with open(args.output, 'w') as f:
        for i in range(len(adjacency_list)):
            f.write(' '.join([str(x) for x in adjacency_list[i]]) + '\n')

