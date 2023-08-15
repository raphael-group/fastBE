import pandas as pd
import numpy as np

import argparse
import json

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
        '-n', '--num-clones', type=int, required=True,
        help='Number of clones in the dataset.'
    )

    parser.add_argument(
        '-o', '--output', type=str, required=True,
        help='Prefix for output files.'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Load with columns as samples and rows as mutations
    variant_matrix = np.loadtxt(args.variant_matrix, dtype=str).T
    total_matrix = np.loadtxt(args.total_matrix, dtype=str).T

    mutation_clustering = pd.read_csv(args.mutation_clustering)

    ssm_file_contents = "id\tname\tvar_reads\ttotal_reads\tvar_read_prob"
    for i in range(variant_matrix.shape[0]):
        out_str = f"s{i}\ts{i}\t"
        out_str += ",".join(variant_matrix[i, :])
        out_str += "\t"
        out_str += ",".join(total_matrix[i, :])
        out_str += "\t"
        # make last column shape[1] 1.0s
        out_str += ",".join(["1.0"] * variant_matrix.shape[1])
        ssm_file_contents += "\n" + out_str 

    with open(f"{args.output}_mutations.ssm", 'w') as f:
        f.write(ssm_file_contents)

    clones_to_mutations = []
    for i in range(args.num_clones):
        muts = mutation_clustering[mutation_clustering['clone'] == i]['mutation'].values.tolist()
        muts = [f"s{j}" for j in muts]
        clones_to_mutations.append(muts)

    samples = [f"sample_{i}" for i in range(variant_matrix.shape[1])]
    with open(f"{args.output}_params.json", 'w') as f:
        data = {
            "samples": samples,
            "clusters": clones_to_mutations,
            "garbage": []
        }
        f.write(json.dumps(data))
