import argparse
import numpy as np
import pandas as pd
import networkx as nx

from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--matrix', type=str, required=False,
                        help='CSV file containing the input matrix B.')
    parser.add_argument('--tree', type=str, required=False,
                        help='Adjacency list describing the input tree.')
    parser.add_argument('--vector', type=str, required=True,
                        help='CSV file containing the input vector f.')

    return parser.parse_args()

def construct_matrix_B(G):
    nodes = list(G.nodes)
    n = len(nodes)

    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if nx.has_path(G, nodes[i], nodes[j]):
                B[j, i] = 1

    return B

def vafpp_linear_program(B, f):
    n = B.shape[0]

    model = Model("1-VAFPP Primal Linear Program")

    u = model.addVars(n, lb=0, vtype=GRB.CONTINUOUS, name="u")
    z = model.addVars(n, lb=0, vtype=GRB.CONTINUOUS, name="z")

    model.setObjective(quicksum(-z[i] for i in range(n)), GRB.MAXIMIZE)

    for i in range(n):
        model.addConstr(z[i] >= f[i] - quicksum(u[j] * B[j, i] for j in range(n)))
        model.addConstr(z[i] >= quicksum(u[j] * B[j, i] for j in range(n)) - f[i])
    model.addConstr(quicksum(u[i] for i in range(n)) <= 1)

    model.optimize()

    if model.status == GRB.OPTIMAL:
        print('Optimal solution found:')
        for i in range(n):
            print('u[{}] = {:.2f}'.format(i, u[i].x))
        for i in range(n):
            print('z[{}] = {:.2f}'.format(i, z[i].x))
    else:
        print('No optimal solution found. Gurobi status code:', model.status)

    u = np.array([u[i].x for i in range(n)])
    print(u.T @ B)
    print(f)
    print(np.sum(np.abs(u.T @ B)))

def solve_dual(tree, f):
    n = len(tree.nodes)

    psi = np.zeros(n)
    gamma = np.zeros(n)

    for i in range(n):
        gamma[i] = max(sum(f[j] for j in tree[i]) - f[i], 0)

    psi = gamma > 0
    obj = np.sum(gamma)

    
    model = Model("1-VAFPP Read Primal From Dual")

    u = model.addVars(n, lb=0, vtype=GRB.CONTINUOUS, name="u")
    z = model.addVars(n, lb=0, vtype=GRB.CONTINUOUS, name="z")

    for i in range(n):
        if psi[i] == 1:
            model.addConstr(u[i] == 0)

        if i == 0:
            if psi[i] == 1:
                model.addConstr(quicksum(u[j] * B[j, i] for j in range(n)) - z[i] == f[i])
            else:
                model.addConstr(z[i] == 0)
                model.addConstr(quicksum(u[j] * B[j, i] for j in range(n)) == f[i])

            continue

        parent = list(tree.predecessors(i))[0]
        
        if (psi[i], psi[parent]) == (1, 1):
            model.addConstr(z[i] == 0)
            model.addConstr(quicksum(u[j] * B[j, i] for j in range(n)) == f[i])
        elif (psi[i], psi[parent]) == (1, 0):
            model.addConstr(quicksum(u[j] * B[j, i] for j in range(n)) - z[i] == f[i])
        elif (psi[i], psi[parent]) == (0, 1):
            model.addConstr(quicksum(u[j] * B[j, i] for j in range(n)) + z[i] == f[i])
        else:
            model.addConstr(z[i] == 0)
            model.addConstr(quicksum(u[j] * B[j, i] for j in range(n)) == f[i])

    model.optimize()

    if model.status == GRB.OPTIMAL:
        print('Optimal solution found:')
        for i in range(n):
            print('u[{}] = {:.2f}'.format(i, u[i].x))
        for i in range(n):
            print('z[{}] = {:.2f}'.format(i, z[i].x))
 
    u = np.array([u[i].x for i in range(n)])
    print(u.T @ B)
    print(f)
    print(np.sum(np.abs(u.T @ B)))

if __name__ == '__main__':
    args = parse_args()

    if args.tree:
        tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
        B = construct_matrix_B(tree)
    else:
        B = pd.read_csv(args.matrix, header=None).values

    f = pd.read_csv(args.vector, header=None).values.flatten()

    vafpp_linear_program(B, f)
    solve_dual(tree, f)
