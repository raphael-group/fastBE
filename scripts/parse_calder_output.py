import numpy as np
import argparse
import networkx as nx
from networkx.drawing.nx_pydot import read_dot

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse CALDER dot file output.')
    parser.add_argument('result', type=str, help='Inferred tree in DOT format.')
    args = parser.parse_args()

    T = read_dot(args.result)

    for node in T.nodes:
        if node == '\\n': continue
        children = ' '.join([n[1:] for n in T[node].keys()])
        print(f'{node[1:]} {children}')
