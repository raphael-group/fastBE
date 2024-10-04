import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Parse the cluster file of Orchard')
    parser.add_argument('clusterings', type=str, help='The input clusterings file')
    parser.add_argument('-k', type=int, required=True, help='The number of clusters')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    clusterings = eval(np.load(args.clusterings)['clusters.json'])
    clustering  = clusterings[-args.k]
    print("mutation,clone")
    for idx, cluster in enumerate(clustering):
        for element in cluster:
            print(f"{element[1:]},{idx}")
