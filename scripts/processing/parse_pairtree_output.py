import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse pairtree output')
    parser.add_argument('result', type=str, help='Result file in NPZ format.')
    parser.add_argument('--output', type=str, help='Output adjacency list.')
    args = parser.parse_args()

    result = np.load(args.result)
    print(list(result['llh'][:10]))
    best_tree_parents = result['struct'][0]

    print(result)
    print(best_tree_parents)
    adjacency_list = [[i] for i in range(len(best_tree_parents) + 1)]
    for i in range(len(best_tree_parents)):
        adjacency_list[best_tree_parents[i]].append(i + 1)
    print(adjacency_list)

    # write adjacency list to NX format
    with open(args.output, 'w') as f:
        for row in adjacency_list[1:]:
            f.write(' '.join([str(x - 1) for x in row]) + '\n')

