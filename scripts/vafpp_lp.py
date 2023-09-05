"""
This file contains several implementations of 1-VAFPP 
solvers. In particular, we implement the following algorithms:
    - A 1-VAFPP solver using linear program
    - A 1-VAFPP solver using the dual linear program
    - A custom 1-VAFPP solver derived via dynamic programming
The first two LP solvers are implemented using Gurobi and are used
for benchmarking the naive algorithm. In contrast, the custom solver
implemented here is a reference implementation for debugging purposes
-- the actual implementation is in C++.
"""

import json
import time
import argparse
import numpy as np
import pandas as pd
import networkx as nx
import piecewise_linear as pl
import matplotlib.pyplot as plt
import cplex

from gurobipy import *

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--clonal-matrix', type=str, required=False,
                        help='Numpy TXT file containing the input matrix B.')
    parser.add_argument('--tree', type=str, required=False,
                        help='Adjacency list describing the input tree.')
    parser.add_argument('--frequency-matrix', type=str, required=True,
                        help='Numpy TXT file containing the input frequency matrix F.')

    parser.add_argument('--output', type=str, required=True,
                        help="Prefix of output files.")

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
def cplex_one_vafpp_linear_program(B, F):
    n, m = F.shape[1], F.shape[0]

    model = cplex.Cplex()
    model.set_problem_name("1-VAFPP Primal Linear Program")

    u_var_indices = list(range(m * n))
    model.variables.add(lb=[0] * (m * n))

    z_var_indices = list(range(m * n, 2 * m * n))
    model.variables.add(lb=[0] * (m * n))

    for i in range(m):
        for j in range(n):
            row_indices = [i * n + k for k in range(n)] + [m * n + i * n + j]
            row_values = [-B[k, j] for k in range(n)] + [1]
            model.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row_indices, val=row_values)], senses="G", rhs=[-F[i, j]])

            row_values = [B[k, j] for k in range(n)] + [1]
            model.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row_indices, val=row_values)], senses="G", rhs=[F[i, j]])

    for i in range(m):
        row_indices = [i * n + j for j in range(n)]
        row_values = [1] * n
        model.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row_indices, val=row_values)], senses="L", rhs=[1])

    obj_indices = list(range(m * n, 2 * m * n))
    model.objective.set_sense(model.objective.sense.minimize)
    model.objective.set_linear(list(zip(obj_indices, [1] * len(obj_indices))))

    starttime = model.get_time()
    model.solve()
    endtime = model.get_time()

    objVal = model.solution.get_objective_value()

    return objVal, endtime - starttime

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

    return model.objVal, model.Runtime

"""
Given a frequency matrix $F$ and a clonal matrix $B$, this function
finds a usage matrix $U$ such that  
        $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$
is minimized, by solving the dual linear program.
"""
def one_vafpp_dual_linear_program(B, F):
    n, m = F.shape[1], F.shape[0]

    model = Model("1-VAFPP Dual Linear Program")
    alpha = model.addMVar(shape=(n, m), lb=0, vtype=GRB.CONTINUOUS, name="alpha")
    beta = model.addMVar(shape=(n, m), lb=0, vtype=GRB.CONTINUOUS, name="beta")
    gamma = model.addMVar(shape=m, lb=0, vtype=GRB.CONTINUOUS, name="gamma")

    model.addConstr(B @ (beta - alpha) + gamma >= 0)
    model.addConstr(alpha + beta <= 1)

    model.setObjective(quicksum(gamma) + quicksum(F[k, i] * (beta[i, k] - alpha[i, k]) for k in range(m) for i in range(n)), GRB.MINIMIZE)
    model.optimize()

    return -1 * model.objVal, model.Runtime

"""
Given a frequency matrix $F$ and a clone tree $T$, this function
finds the minimizing value of
        $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$
over all usage matrices in $\mathcal{O}(mn^2)$ using dynamic
programming and an analysis of convex and continuous, 
piecewise linear functions.
"""
def one_vafpp_dp(tree, F):
    W = np.zeros(F.shape)
    for j in range(F.shape[0]):
        for i in range(F.shape[1]):
            W[j, i] = F[j, i] - sum(F[j, k] for k in tree[i]) 

    def one_vafpp_dual_recursive(j, i):
        if len(tree[i]) == 0:
            return pl.PiecewiseLinear([W[j, i]], 0) 

        g_out = pl.PiecewiseLinear([W[j, i]], 0)
        for k in tree[i]:
            f = one_vafpp_dual_recursive(j, k)
            g = pl.compute_minimizer(f)
            g_out = g_out + g

        return g_out

    obj = 0
    for j in range(F.shape[0]):
        f = one_vafpp_dual_recursive(j, 0)
        f = pl.compute_minimizer(f) + pl.PiecewiseLinear([1 - F[j, 0]], 0)
        row_obj = np.min(f.intercepts)
        obj += row_obj
        
    return -1 * obj

if __name__ == '__main__':
    args = parse_args()

    if args.tree:
        tree = nx.read_adjlist(args.tree, nodetype=int, create_using=nx.DiGraph())
        B = construct_clonal_matrix(tree)
    else:
        B = np.loadtxt(args.clonal_matrix)

    F = np.loadtxt(args.frequency_matrix)

    results = {}

    start = time.time()
    obj1, gurobi_time = one_vafpp_linear_program(B, F)
    end = time.time()

    results["gurobi_lp_time_with_building"] = end - start
    results["gurobi_lp_time_without_building"] = gurobi_time
    results["gurobi_lp_obj"] = obj1

    start = time.time()
    obj2, gurobi_time = one_vafpp_dual_linear_program(B, F)
    end = time.time()

    results["gurobi_dual_time_with_building"] = end - start
    results["gurobi_dual_time_without_building"] = gurobi_time
    results["gurobi_dual_obj"] = obj2

    start = time.time()
    obj3 = one_vafpp_dp(tree, F)
    end = time.time()

    results["dp_time"] = end - start
    results["dp_obj"] = obj3

    start = time.time()
    obj4, cplex_time = cplex_one_vafpp_linear_program(B, F)
    end = time.time()

    results["cplex_lp_time_with_building"] = end - start
    results["cplex_lp_time_without_building"] = cplex_time
    results["cplex_lp_obj"] = obj4

    with open(f"{args.output}_results.json", 'w') as f:
        json.dump(results, f, indent=4)
