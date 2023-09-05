import argparse
import json
import os
import pandas as pd

def read_and_process(directory, res_type):
    data = []
    for filename in os.listdir(directory):
        if not res_type in filename:
            continue

        if filename.endswith('.json'):
            with open(os.path.join(directory, filename), 'r') as f:
                file_data = json.load(f)
                if 'dual_obj' in file_data:
                    del file_data['dual_obj']  # exclude 'dual_obj'
                parameters = parse_filename(filename, res_type)
                data.append({**parameters, **file_data})

    df = pd.DataFrame(data)
    return df

def parse_filename(filename, res_type):
    base_name = filename.split(f'_{res_type}_results.json')[0]
    params = base_name.split('_')
    parameters = {
        'm': int(params[0][1:]), 
        'n': int(params[1][1:]), 
        's': int(params[2][1:]), 
        'c': int(params[3][1:]), 
        'r': int(params[4][1:])
    }
    return parameters

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some json files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory where the json files are located.')
    parser.add_argument('--output', type=str, help='The output file name.')
    parser.add_argument('--type', choices=['python', 'cpp'], default='python')
    
    args = parser.parse_args()

    df = read_and_process(args.directory, args.type)
    if args.output:
        df.to_csv(args.output, index=False)
