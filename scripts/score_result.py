import json
import argparse
import networkx as nx
import numpy as np

from gurobipy import *

"""
Given a frequency matrix $F$ and a clonal matrix $B$, this function
finds a usage matrix $U$ such that  
        $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$
is minimized, by using linear programming.
"""
def one_vafpp_linear_program(B, F):
    n, m = F.shape[1], F.shape[0]

    model = Model("1-VAFPP Primal Linear Program")

    U = model.addMVar(shape=(m, n), lb=0, vtype=GRB.CONTINUOUS, name="U")
    Z = model.addMVar(shape=(m, n), lb=0, vtype=GRB.CONTINUOUS, name="Z")

    model.addConstr(Z >= F - U @ B)
    model.addConstr(Z >= U @ B - F)

    for k in range(m):
        model.addConstr(U[k, :].sum() <= 1)

    model.setObjective(Z.sum(), GRB.MINIMIZE)
    model.optimize()

    # convert U to a numpy array
    U_np = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            U_np[i, j] = U[i, j].x

    return U_np

"""
Constructs the clonal matrix from a clonal tree.
"""
def construct_clonal_matrix(tree):
    n = len(tree.nodes)

    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if nx.has_path(tree, i, j):
                B[j, i] = 1

    return B

"""
Computes the set of all pairwise relations in a directed tree.
"""
def get_relations(tree):
    relations = set()
    for u in tree.nodes():
        for v in tree.nodes():
            if u == v:
                continue

            if nx.has_path(tree, u, v):
                relations.add((u, v))

    return relations

def parse_args():
    parser = argparse.ArgumentParser(description="Score ")
    parser.add_argument('true_tree', help='Adjacency list of true tree')
    parser.add_argument('true_usage_matrix', help='Frequency matrix of true tree')
    parser.add_argument('inferred_tree', help='Adjacency list of inferred tree')
    parser.add_argument('-o', '--output', help='Name of JSON output file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    true_tree = nx.read_adjlist(args.true_tree, nodetype=int, create_using=nx.DiGraph)
    inferred_tree = nx.read_adjlist(args.inferred_tree, nodetype=int, create_using=nx.DiGraph)

    # Compute the pairwise relation error
    true_relations = get_relations(true_tree)
    inferred_relations = get_relations(inferred_tree)

    false_positives = inferred_relations - true_relations
    positives = true_relations
    false_negatives = true_relations - inferred_relations
    negatives = set([(u, v) for u in true_tree.nodes() for v in true_tree.nodes() if u != v]) - positives

    fpr = len(false_positives) / len(positives)

    if len(negatives) == 0:
        fnr = 0
    else:
        fnr = len(false_negatives) / len(negatives)

    result = {
        'pairwise_relations': {
            'false_positive_rate': fpr,
            'false_negative_rate': fnr,
            'false_positives': len(false_positives),
            'false_negatives': len(false_negatives),
            'positives': len(positives),
            'negatives': len(negatives)
        }
    }

    # Compute the frequency and usage matrix error
    U = np.loadtxt(args.true_usage_matrix)
    B = construct_clonal_matrix(true_tree)
    F = U @ B

    B_hat = construct_clonal_matrix(inferred_tree)
    U_hat = one_vafpp_linear_program(B_hat, F)

    u_l1_error = np.sum(np.abs(U - U_hat))
    f_l1_error = np.sum(np.abs(F - U @ B_hat))

    result['frequency_matrix'] = {
        'F_error': f_l1_error,
        'U_error': u_l1_error
    }

    # Write the result to a JSON file 
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f)
    else:
        print(json.dumps(result, indent=4))

