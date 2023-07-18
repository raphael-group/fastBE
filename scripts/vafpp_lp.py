import argparse
import numpy as np
import pandas as pd
import networkx as nx

from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--clonal-matrix', type=str, required=False,
                        help='Numpy TXT file containing the input matrix B.')
    parser.add_argument('--tree', type=str, required=False,
                        help='Adjacency list describing the input tree.')
    parser.add_argument('--frequency-matrix', type=str, required=True,
                        help='Numpy TXT file containing the input frequency matrix F.')

    return parser.parse_args()

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

    for i in range(n):
        for k in range(m):
            model.addConstr(Z[k, i] >= F[k, i] - U[k, :] @ B[:, i])
            model.addConstr(Z[k, i] >= U[k, :] @ B[:, i] - F[k, i])

    for k in range(m):
        model.addConstr(U[k, :].sum() <= 1)

    model.setObjective(Z.sum(), GRB.MINIMIZE)
    model.optimize()

    # print U as a matrix
    print("U:")
    for i in range(m):
        for j in range(n):
            print(U[i, j].x, end=" ")
        print()

    return model.objVal

"""
Given a frequency matrix $F$ and a clone tree $T$, this function
finds the minimizing value of
        $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$
over all usage matrices.
"""
def one_vafpp_dual(tree, F):
    n, m = F.shape[1], F.shape[0]

    psi = np.zeros((m, n))
    gamma = np.zeros((m, n))

    for i in range(n):
        for k in range(m):
            gamma[k, i] = max(sum(F[k, j] for j in tree[i]) - F[k, i], 0)

    psi = gamma > 0
    obj = np.sum(gamma)
    return obj

if __name__ == '__main__':
    args = parse_args()

    if args.tree:
        tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
        B = construct_clonal_matrix(tree)
    else:
        B = np.loadtxt(args.clonal_matrix)

    F = np.loadtxt(args.frequency_matrix)

    obj1 = one_vafpp_linear_program(B, F)
    obj2 = one_vafpp_dual(tree, F)

    print("Linear Program Objective: {}".format(obj1))
    print("Dual Objective: {}".format(obj2))

