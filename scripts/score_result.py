import json
import argparse
import networkx as nx

def parse_args():
    parser = argparse.ArgumentParser(description="Score ")
    parser.add_argument('true_tree', help='Adjacency list of true tree')
    parser.add_argument('inferred_tree', help='Adjacency list of inferred tree')
    return parser.parse_args()

def get_relations(tree):
    relations = set()
    for u in tree.nodes():
        for v in tree.nodes():
            if u == v:
                continue

            if nx.has_path(tree, u, v):
                relations.add((u, v))

    return relations

if __name__ == "__main__":
    args = parse_args()

    true_tree = nx.read_adjlist(args.true_tree, create_using=nx.DiGraph)
    inferred_tree = nx.read_adjlist(args.inferred_tree, create_using=nx.DiGraph)
        
    true_relations = get_relations(true_tree)
    inferred_relations = get_relations(inferred_tree)

    # set diff
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

    print(json.dumps(result))
