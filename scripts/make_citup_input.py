import sys
import argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'frequency_matrix', type=str,
        help='Frequency matrix in TXT format.'
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    frequency_matrix = np.loadtxt(args.frequency_matrix).T
    np.savetxt(sys.stdout, frequency_matrix, fmt='%.6f')
