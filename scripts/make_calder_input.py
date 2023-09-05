import pandas as pd
import numpy as np

import argparse
import json

def flatten(l):
    return [item for sublist in l for item in sublist]

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'variant_matrix', type=str,
        help='Variant read count matrix in TXT format.'
    )

    parser.add_argument(
        'total_matrix', type=str,
        help='Total read count matrix in TXT format.'
    )

    parser.add_argument(
        'mutation_clustering', type=str,
        help='Mapping from mutations to clones in CSV format.'
    )

    parser.add_argument(
        '-o', '--output', type=str, required=True,
        help='Prefix for output files.'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Load with columns as samples and rows as mutations
    variant_matrix = np.loadtxt(args.variant_matrix, dtype=int).T
    total_matrix = np.loadtxt(args.total_matrix, dtype=int).T

    mutation_clustering = pd.read_csv(args.mutation_clustering)
    num_clones = mutation_clustering['clone'].max() + 1

    clones_to_mutations = []
    for i in range(num_clones):
        muts = mutation_clustering[mutation_clustering['clone'] == i]['mutation'].values.tolist()
        clones_to_mutations.append(muts)

    columns = []
    for i in range(num_clones):
        muts = clones_to_mutations[i]
        columns.append(total_matrix[muts].sum(axis=0))
        columns.append(variant_matrix[muts].sum(axis=0))

    column_names = flatten([[f'm{i}', f'm{i}'] for i in range(num_clones)])
    df = pd.DataFrame(columns).T
    df.columns = column_names
    df.to_csv(args.output, index=True, sep='\t')
