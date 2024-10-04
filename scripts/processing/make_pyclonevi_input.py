import sys
import numpy as np
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make PyCloneVI input from variant and total read matrices')
    parser.add_argument('-v', '--variant_matrix', type=str, help='Path to variant matrix')
    parser.add_argument('-t', '--total_read_matrix', type=str, help='Path to total read matrix')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    variant_matrix = np.loadtxt(args.variant_matrix)
    total_read_matrix = np.loadtxt(args.total_read_matrix)

    rows = []
    m = variant_matrix.shape[0]
    n = variant_matrix.shape[1]
    for i, j in np.ndindex(m, n):
        rows.append({
            'ref_counts': total_read_matrix[i, j] - variant_matrix[i, j],
            'alt_counts': variant_matrix[i, j],
            'mutation_id': 'm' + str(j),
            'sample_id': 's' + str(i),
            'major_cn': 2,
            'minor_cn': 2,
            'normal_cn': 2,
        })

    df = pd.DataFrame(rows)
    df.to_csv(sys.stdout, sep='\t', index=False)
