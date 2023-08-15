import argparse
import numpy as np
import pandas as pd
import networkx as nx

from scipy.stats import poisson, binom

"""
Simulate a clonal tree with n nodes and m mutations,
by assigning each of the m mutations to a random node.

Output:
    - tree: a networkx DiGraph object representing the clonal tree
      where the root is labeled by 0, and the other nodes are labeled
      by 1, 2, ..., n - 1. The attached mutations are stored in the 
      'mutation' attribute of the nodes.
"""
def simulate_clonal_tree(m, n):
    assert m >= n

    tree = nx.DiGraph()
    tree.add_node(0) 

    for i in range(1, n):
        parent = np.random.choice(np.arange(i))
        tree.add_node(i)
        tree.add_edge(parent, i)

    # Ensures that every node has at least one mutation
    remaining_nodes = list(range(n))
    for i in range(n):
        node = np.random.choice(list(remaining_nodes))
        remaining_nodes.remove(node)

        if 'mutation' not in tree.nodes[node]:
            tree.nodes[node]['mutation'] = []

        tree.nodes[node]['mutation'].append(i)

    for i in range(n, m):
        node = np.random.choice(np.arange(n))

        if 'mutation' not in tree.nodes[node]:
            tree.nodes[node]['mutation'] = []

        tree.nodes[node]['mutation'].append(i)

    mutation_to_clone_mapping = {}
    for node in tree.nodes:
        if 'mutation' not in tree.nodes[node]:
            continue

        for mutation in tree.nodes[node]['mutation']:
            mutation_to_clone_mapping[mutation] = node

    return tree, mutation_to_clone_mapping

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
Simulate a usage matrix with s samples and n clones,
by randomly selecting k of n clones and assigning each of them
a random usage probability.

Output:
    - matrix: a s x n matrix where each row represents a sample 
    and the sum of the entries is 1.
"""
def simulate_usage_matrix(tree, s, n):
    matrix = np.zeros((s, len(tree)))

    for i in range(s):
        clones = np.random.choice(np.arange(len(tree)), size=np.random.randint(1, n), replace=False)
        usages = np.random.dirichlet(np.ones(len(clones)), size=1)[0]
        for clone, usage in zip(clones, usages):
            matrix[i, clone] = usage

    return matrix

"""
Simulate a mutation read counts from a clonal matrix, usage matrix and 
a mapping from mutations to clone.
"""
def simulate_read_counts(usage_matrix, clonal_matrix, mutation_to_clone_mapping, num_mutations, coverage):
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
                    F[s, mutation_to_clone_mapping[mutation]]
            )
            # there could be a bug here...
            variant_count_matrix[s, mutation] = f

    assert np.all(variant_count_matrix <= total_count_matrix)
    return variant_count_matrix, total_count_matrix

def observe_frequency_matrix(variant_count_matrix, total_count_matrix, mutation_to_clone_mapping, num_clones):
    clone_to_mutation_mapping = {}
    for clone in range(num_clones):
        clone_to_mutation_mapping[clone] = []

    for mutation in mutation_to_clone_mapping:
        clone = mutation_to_clone_mapping[mutation]
        clone_to_mutation_mapping[clone].append(mutation)

    clone_mut = lambda c: clone_to_mutation_mapping[c]

    obs_frequency_matrix = np.zeros((variant_count_matrix.shape[0], len(clone_to_mutation_mapping.keys())))
    for s in range(obs_frequency_matrix.shape[0]):
        for clone in range(obs_frequency_matrix.shape[1]):
            variant_reads = sum([variant_count_matrix[s, m] for m in clone_mut(clone)])
            total_reads   = sum([total_count_matrix[s, m] for m in clone_mut(clone)])
            if total_reads > 0:
                obs_frequency_matrix[s, clone] = variant_reads / total_reads

    return obs_frequency_matrix

def main():
    parser = argparse.ArgumentParser(description='Simulate a clonal matrix, usage matrix and read counts.')

    parser.add_argument('--mutations', type=int, required=True, help='Number of mutations.')
    parser.add_argument('--clones', type=int, required=True, help='Number of clones.')
    parser.add_argument('--samples', type=int, required=True, help='Number of sequenced samples.')
    parser.add_argument('--coverage', type=float, required=True, help='Expected sequencing coverage.')
    parser.add_argument('--output', type=str, required=True, help='Output prefix.')
    parser.add_argument('--seed', type=int, default=0, help='Random seed.')

    args = parser.parse_args()

    assert args.mutations >= args.clones, 'Number of mutations must be greater than or equal to number of clones.'

    np.random.seed(args.seed)

    tree, mutation_to_clone_mapping = simulate_clonal_tree(args.mutations, args.clones)
    clonal_matrix = construct_clonal_matrix(tree)
    usage_matrix = simulate_usage_matrix(tree, args.samples, args.clones)
    variant_matrix, total_matrix = simulate_read_counts(
            usage_matrix, clonal_matrix, mutation_to_clone_mapping, 
            args.mutations, args.coverage
    )

    f_hat = observe_frequency_matrix(variant_matrix, total_matrix, mutation_to_clone_mapping, args.clones)
    
    np.savetxt(f'{args.output}_clonal_matrix.txt', clonal_matrix, fmt='%d')
    np.savetxt(f'{args.output}_usage_matrix.txt', usage_matrix, fmt='%.4f')
    np.savetxt(f'{args.output}_variant_matrix.txt', variant_matrix, fmt='%d')
    np.savetxt(f'{args.output}_total_matrix.txt', total_matrix, fmt='%d')
    np.savetxt(f'{args.output}_obs_frequency_matrix.txt', f_hat, fmt='%.4f')

    nx.write_adjlist(tree, f'{args.output}_tree.txt')
    
    df = pd.DataFrame.from_dict(mutation_to_clone_mapping, orient='index').reset_index()
    df.columns = ['mutation', 'clone']
    df.to_csv(f'{args.output}_mutation_to_clone_mapping.txt', sep=',', index=False)

if __name__ == "__main__":
    main()

