import numpy as np
import pandas as pd
import argparse
import json 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('somatic_mutations', help='Somatic mutations (.ssm) file', type=str)
    parser.add_argument('params_json', help='JSON file containing parameters for Pairtree.', type=str)
    parser.add_argument('--output', help='Prefix for output files.', type=str, required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    somatic_mutations = pd.read_csv(args.somatic_mutations, sep='\t')
    with open(args.params_json, 'r') as f:
        params = json.load(f)

    somatic_mutations = somatic_mutations.set_index('id')
    mutation_ids = somatic_mutations.index.tolist()

    variant_matrix_rows = []
    total_matrix_rows = []
    for mutation_id in mutation_ids:
        variant_matrix_row = somatic_mutations.loc[mutation_id].var_reads.split(',')
        variant_matrix_row = np.array(list(map(int, variant_matrix_row)))
        total_matrix_row = somatic_mutations.loc[mutation_id].total_reads.split(',')
        total_matrix_row = np.array(list(map(int, total_matrix_row)))
        variant_matrix_rows.append(variant_matrix_row)
        total_matrix_rows.append(total_matrix_row)

    variant_matrix = np.vstack(variant_matrix_rows).T
    total_matrix = np.vstack(total_matrix_rows).T

    np.savetxt(f'{args.output}_variant_matrix.txt', variant_matrix, fmt='%d')
    np.savetxt(f'{args.output}_total_matrix.txt', total_matrix, fmt='%d')
    np.savetxt(f'{args.output}_samples.txt', params['samples'], fmt='%s')
    np.savetxt(f'{args.output}_mutations.txt', mutation_ids, fmt='%s')

    mutation_idx = {j : i for i, j in enumerate(mutation_ids)}
    mutation_clone_mapping = {}
    for (i, cluster) in enumerate(params['clusters']):
        for mutation_id in cluster:
            mutation_clone_mapping[mutation_idx[mutation_id]] = i

    clone_mutation_mapping = {}
    for mutation in mutation_clone_mapping.keys():
        clone = mutation_clone_mapping[mutation]
        if clone not in clone_mutation_mapping:
            clone_mutation_mapping[clone] = []
        clone_mutation_mapping[clone].append(mutation)

    mutation_clone_df = pd.DataFrame.from_dict(mutation_clone_mapping, orient='index').reset_index()
    mutation_clone_df.columns = ['mutation', 'clone']
    mutation_clone_df.to_csv(f'{args.output}_mutation_clone_mapping.csv', index=False)

    frequency_matrix_cols = []
    for clone in clone_mutation_mapping.keys():
        variant_col_sum = np.zeros(variant_matrix.shape[0])
        total_col_sum = np.zeros(total_matrix.shape[0])
        for mutation in clone_mutation_mapping[clone]:
            variant_col_sum += variant_matrix[:, mutation]
            total_col_sum += total_matrix[:, mutation]

        frequency_matrix_cols.append(variant_col_sum / total_col_sum)

    frequency_matrix = np.vstack(frequency_matrix_cols).T

    np.savetxt(f'{args.output}_frequency_matrix.txt', frequency_matrix, fmt='%.4f')
    print(frequency_matrix.shape)
