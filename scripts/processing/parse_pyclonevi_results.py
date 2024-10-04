import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Parse PyCloneVI results')
    parser.add_argument('results', type=str, help='PyCloneVI results file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    results = pd.read_csv(args.results, sep='\t')

    print("mutation,clone")
    for idx, mutations in enumerate(results.groupby('cluster_id').mutation_id.unique()):
        for mutation in mutations:
            print(f"{mutation[1:]},{idx}")

