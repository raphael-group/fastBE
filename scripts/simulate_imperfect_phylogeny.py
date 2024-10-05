#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import networkx as nx

from scipy.stats import poisson, binom

def rand_int(rng, a, b):
    return rng.integers(a, b+1)

def rand_predecessor(node, predecessors, weights, rng):
    weighted_preds = [(p, weights[(p, node)]) for p in predecessors]
    total_weight = sum(weight for _, weight in weighted_preds)
    r = rng.random() * total_weight
    upto = 0
    for pred, weight in weighted_preds:
        if upto + weight >= r:
            return pred
        upto += weight

def sample_random_spanning_tree(G, weights, rng, root=None):
    spanning_tree = nx.DiGraph()

    for u in G.nodes:
        spanning_tree.add_node(u)

    next_node = [-1] * len(G.nodes)
    in_tree = [False] * len(G.nodes)

    if root is None:
        root = rand_int(rng, 0, len(G.nodes) - 1)

    in_tree[root] = True
    for u in G.nodes:
        if in_tree[u]:
            continue

        v = u
        while not in_tree[v]:
            pred = list(G.predecessors(v))
            if len(pred) == 0:
                raise RuntimeError("Graph is not strongly connected")

            next_node[v] = rand_predecessor(v, pred, weights, rng)
            v = next_node[v]

        v = u
        while not in_tree[v]:
            in_tree[v] = True
            spanning_tree.add_edge(next_node[v], v)
            v = next_node[v]

    return spanning_tree, root

"""
Simulate a clonal tree with n mutations.

Output:
    - tree: a networkx DiGraph object representing the clonal tree
      where the root is labeled by 0, and the other nodes are labeled
      by 1, 2, ..., n - 1. 
"""
def simulate_clonal_tree(n, seed):
    # construct complete graph
    G = nx.DiGraph()
    for i in range(n):
        G.add_node(i)

    for i in range(n):
        for j in range(n):
            if i != j:
                G.add_edge(i, j)

    # sample random spanning tree
    rng = np.random.default_rng(seed)
    tree, _ = sample_random_spanning_tree(G, {(i, j): 1 for i in range(n) for j in range(n) if i != j}, rng, root=0)

    return tree

"""
Constructs the clonal matrix from a clonal tree.
"""
def construct_clonal_matrix(tree, mutations, node_to_mutation):
    n = len(tree.nodes)

    rows = []
    for i in range(n):
        ancestors = list(nx.shortest_path(tree, 0, i))
        muts = [node_to_mutation[node] for node in ancestors]
        row = np.zeros(len(mutations))
        for mut in muts:
            row[mut[0]] = 1 if mut[1] == '+' else 0
        rows.append(row)
        
    B = np.array(rows)
    return B

"""
Simulate a usage matrix with s samples and n clones,
by randomly selecting k of n clones and assigning each of them
a random usage probability.

Output:
    - matrix: a s x n matrix where each row represents a sample 
    and the sum of the entries is 1.
"""
def simulate_usage_matrix(tree, s):
    n = len(tree)
    matrix = np.zeros((s, len(tree)))

    for i in range(s):
        clones = np.random.choice(np.arange(len(tree)), size=np.random.randint(1, n), replace=False)
        usages = np.random.dirichlet(np.ones(len(clones)), size=1)[0]
        for clone, usage in zip(clones, usages):
            matrix[i, clone] = usage

    return matrix

"""
# TODO: update this part of the code
Simulate a mutation read counts from a clonal matrix, usage matrix and 
a mapping from mutations to clone.
"""
def simulate_read_counts(usage_matrix, clonal_matrix, num_mutations, coverage):
    F = usage_matrix @ clonal_matrix

    # Clamp to [0, 1] for floating point errors
    F[F < 0] = 0
    F[F > 1] = 1
    
    variant_count_matrix = np.zeros((usage_matrix.shape[0], num_mutations))
    total_count_matrix   = np.zeros((usage_matrix.shape[0], num_mutations))
    total_count_matrix   = poisson.rvs(coverage, size=total_count_matrix.shape)
    for mutation in range(num_mutations):
        for s in range(usage_matrix.shape[0]):
            f = binom.rvs(
                    total_count_matrix[s, mutation],
                    F[s, mutation]
            )
            # there could be a bug here...
            variant_count_matrix[s, mutation] = f

    assert np.all(variant_count_matrix <= total_count_matrix)
    return variant_count_matrix, total_count_matrix

def main():
    parser = argparse.ArgumentParser(description='Simulate a clonal matrix, usage matrix and read counts.')

    parser.add_argument('--mutations', type=int, required=True, help='Number of mutations.')
    parser.add_argument('--violations', type=int, default=0, help='Number of violations.')
    parser.add_argument('--samples', type=int, required=True, help='Number of sequenced samples.')
    parser.add_argument('--coverage', type=float, required=True, help='Expected sequencing coverage.')
    parser.add_argument('--output', type=str, required=True, help='Output prefix.')
    parser.add_argument('--seed', type=int, default=0, help='Random seed.')

    args = parser.parse_args()

    np.random.seed(args.seed)

    tree = simulate_clonal_tree(args.mutations + args.violations, args.seed)
    mutations = list(range(args.mutations))

    # classifies mutations as gain or loss
    node_to_mutation = {}
    for i in range(args.mutations):
        node_to_mutation[i] = (i, '+')

    # randomly assign gain or loss to violations
    for i in range(args.violations):
        node = args.mutations + i
        ancestors = list(nx.ancestors(tree, node))

        if len(ancestors) <= 1:
            mut = np.random.choice(list(set(mutations) - set([0])))
            node_to_mutation[node] = (mut, '+')
            continue

        loss = np.random.choice([True, False])
        if loss:
            mut = np.random.choice(ancestors)
            node_to_mutation[node] = (mut, '-')
        else:
            mut = np.random.choice(list(set(mutations) - set(ancestors)))
            node_to_mutation[node] = (mut, '+')

    clonal_matrix = construct_clonal_matrix(tree, mutations, node_to_mutation)
    usage_matrix = simulate_usage_matrix(tree, args.samples)
    variant_matrix, total_matrix = simulate_read_counts(
            usage_matrix, clonal_matrix,
            args.mutations, args.coverage
    )

    f_hat = variant_matrix / total_matrix
    
    node_to_mutation = pd.DataFrame(node_to_mutation, index=['mutation', 'type']).T.reset_index()
    node_to_mutation.columns = ['node', 'mutation', 'type']
    np.savetxt(f'{args.output}_clonal_matrix.txt', clonal_matrix, fmt='%d')
    np.savetxt(f'{args.output}_usage_matrix.txt', usage_matrix, fmt='%.4f')
    np.savetxt(f'{args.output}_frequency_matrix.txt', f_hat, fmt='%.4f')
    np.savetxt(f'{args.output}_variant_matrix.txt', variant_matrix, fmt='%d')
    np.savetxt(f'{args.output}_total_matrix.txt', total_matrix, fmt='%d')
    node_to_mutation.to_csv(f'{args.output}_node_to_mutation.txt', sep='\t', index=False)
    nx.write_adjlist(tree, f'{args.output}_tree.txt')
    
if __name__ == "__main__":
    main()
