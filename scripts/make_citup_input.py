import sys
import argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'frequency_matrix', type=str,
        help='Frequency matrix in TXT format.'
    )

    parser.add_argument(
        '--output', type=str, default='citup_'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    frequency_matrix = np.loadtxt(args.frequency_matrix).T
    np.savetxt(args.output + 'frequency_matrix.txt', frequency_matrix, fmt='%.6f')
    with open(args.output + 'clustering.txt', 'w') as f:
        for i in range(frequency_matrix.shape[0]):
            f.write(str(i) + '\n')
