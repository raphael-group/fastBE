import re
import os
import json
import argparse
import pandas as pd

def load_files(directory):
    data = []
    for file in os.listdir(directory):
        if file.endswith(".json"):
            with open(os.path.join(directory, file), 'r') as f:
                content = json.load(f)
                match = re.search(r'r(\d+)_m(\d+)_n(\d+)_s(\d+)_c(\d+)_r(\d+)', file)
                iteration, m, n, s, c, r = match.groups()

                timing_result_file = os.path.join(f'nextflow_results/benchmarking/allele_minima/', f'r{iteration}_m{m}_n{n}_s{s}_c{c}_r{r}_timing.txt')
                if os.path.exists(timing_result_file):
                    with open(timing_result_file, 'r') as f:
                        timing = f.read()
                        match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)', timing) 
                        elapsed_time = match.groups()[0]
                        elapsed_time = sum(x * float(t) for x, t in zip([1, 60, 3600], elapsed_time.split(":")[::-1]))
                else:
                    elapsed_time = None

                row = {
                    'iteration': iteration,
                    'mutations': m,
                    'clones': n,
                    'samples': s,
                    'coverage': c,
                    'seed': r,
                    'objective_value': content['objective_value'],
                    'elapsed_time': elapsed_time
                }
                data.append(row)

    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Load JSON files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory containing the JSON files.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('benchmark_summary.csv', index=False)

if __name__ == '__main__':
    main()

