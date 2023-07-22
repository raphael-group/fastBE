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

    print(B)
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
            model.addConstr(Z[k, i] >= F[k, i] - quicksum(U[k, j] * B[j, i] for j in range(n)))
            model.addConstr(Z[k, i] >= quicksum(U[k, j] * B[j, i] for j in range(n)) - F[k, i])

    for k in range(m):
        model.addConstr(U[k, :].sum() <= 1)

    model.setObjective(Z.sum(), GRB.MINIMIZE)
    model.optimize()

    # print model constraints without commenting out lines:
    print(model.display())

    U_m = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            U_m[i, j] = U[i, j].x

    print(f"U = {U_m}")
    print(f"UB = {U_m @ B}")
    print(f"F - UB = {F - U_m @ B}")
    print(f"F = {F}")

    return model.objVal

# DUAL APPEARS TO WORK
def one_vafpp_dual_linear_program(B, F):
    n = F.shape[1]
    F = F[0, :]

    model = Model("1-VAFPP Dual Linear Program")

    alpha = model.addMVar(shape=n, lb=0, vtype=GRB.CONTINUOUS, name="alpha")
    beta = model.addMVar(shape=n, lb=0, vtype=GRB.CONTINUOUS, name="beta")
    gamma = model.addMVar(shape=1, lb=0, vtype=GRB.CONTINUOUS, name="gamma")

    for i in range(n):
        model.addConstr(quicksum(B[i, j] * (beta[j] - alpha[j]) for j in range(n)) + gamma[0] >= 0)
        model.addConstr(alpha[i] + beta[i] <= 1)

    model.setObjective(gamma[0] + quicksum(F[i] * (beta[i] - alpha[i]) for i in range(n)), GRB.MINIMIZE)
    model.optimize()

    # print model constraints without commenting out lines:
    print(model.display())

def one_vafpp_dual_linear_program2(B, F):
    n = F.shape[1]

    model = Model("1-VAFPP Dual Linear Program")

    lambda_var = model.addMVar(shape=n, lb=-1, ub=1, vtype=GRB.CONTINUOUS, name="lambda")
    gamma = model.addMVar(shape=1, lb=0, vtype=GRB.CONTINUOUS, name="gamma")
    psi = model.addMVar(shape=n, lb=0, vtype=GRB.CONTINUOUS, name="psi")

    model.addConstr(B @ lambda_var == psi - gamma[0] * np.ones(n))

    model.setObjective(gamma[0] + F @ lambda_var, GRB.MINIMIZE)
    model.optimize()

    # print model constraints without commenting out lines:
    print(model.display())
    return model.objVal

def one_vafpp_dual_linear_program3(f, tree): 
    n = len(f)

    model = Model("1-VAFPP Linear Program")

    psi = model.addMVar(shape=n, lb=0, vtype=GRB.CONTINUOUS, name="psi")
    gamma = model.addMVar(shape=1, lb=0, vtype=GRB.CONTINUOUS, name="gamma")

    objective = gamma[0] * (1 - f[0])

    for i in range(n):
        child_sum = sum(f[j] for j in tree[i])
        objective += psi[i] * (f[i] - child_sum)

        if i != 0:
            parent_j = list(tree.predecessors(i))[0]
            model.addConstr(psi[i] - psi[parent_j] <= 1)
            model.addConstr(psi[i] - psi[parent_j] >= -1)
        else:
            model.addConstr(psi[i] - gamma[0] <= 1)
            model.addConstr(psi[i] - gamma[0] >= -1)

    model.setObjective(objective, GRB.MINIMIZE)
    model.optimize()

    print(model.display())
    print("psi = {}".format(psi.x))
    print("gamma = {}".format(gamma.x))

    return model.objVal

def one_vafpp_dual_recursive(tree, w, i):
    if len(tree[i]) == 0:
        print(w[i])
        return (w[i], 0) # returns (a, b) as defined in algo

    children = [one_vafpp_dual_recursive(tree, w, j) for j in tree[i]]
    a_out = w[i] + sum(a for (a, _) in children)
    print("i = {}, tree[i] = {}".format(i, [j for j in tree[i]]))
    print("a_out = {}".format(a_out))
    print(children)
    b_out = sum(b - np.abs(a) for (a, b) in children)
    return (a_out, b_out)

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
            print("i = {}, tree[i] = {}".format(i, [j for j in tree[i]]))
            gamma[k, i] = max(sum(F[k, j] for j in tree[i]) - F[k, i], 0)
            print("gamma[{}, {}] = {}".format(k, i, gamma[k, i]))

    psi = gamma > 0
    obj = np.sum(gamma)

    print(obj)
    f = F[0, :]
    psi = psi[0, :]
    gamma = gamma[0, :]

    print(psi)
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
    return obj

if __name__ == '__main__':
    args = parse_args()

    if args.tree:
        tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
        B = construct_clonal_matrix(tree)
    else:
        B = np.loadtxt(args.clonal_matrix)

    F = np.loadtxt(args.frequency_matrix)[[13], :]

    print(F)

    obj1 = one_vafpp_linear_program(B, F)
    obj2 = one_vafpp_dual_linear_program2(B, F)
    obj2 = one_vafpp_dual_linear_program3(F[0, :], tree)

    w = np.zeros(F.shape[1])
    for i in range(F.shape[1]):
        w[i] = F[0, i] - sum(F[0, j] for j in tree[i]) 

    obj4 = one_vafpp_dual_recursive(tree, w, 0)
    print(obj4)

    print("Linear Program Objective: {}".format(obj1))
    print("Dual Objective: {}".format(obj2))

