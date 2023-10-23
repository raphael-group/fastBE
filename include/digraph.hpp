#ifndef _DIGRAPH_H
#define _DIGRAPH_H

#include <map>
#include <unordered_map>
#include <random>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

/*
  Simple adjacency list representation of a DiGraph where the 
  vertices are represented as integers from 0 to N - 1. And
  the vertices have data associated with them.
*/

template <class T>
class vertex {
public:
    int id;
    T data;
    vertex() {};
    vertex(int id, T data) : id(id), data(data) {};
};

template <class T>
class digraph {
private:
    int id_counter = 0;

    std::unordered_map<int, std::set<int>> succ;
    std::unordered_map<int, std::set<int>> pred;

    std::unordered_map<int, vertex<T>> vertices;
public:
    // returns id of created vertex
    int add_vertex(T data) {
        vertex<T> v(id_counter, data);
        vertices[v.id] = v;
        succ[v.id] = std::set<int>();
        pred[v.id] = std::set<int>();
        id_counter++;
        return v.id;
    }

    void add_edge(int u, int v) {
        succ[u].insert(v);
        pred[v].insert(u);
    }

    void remove_edge(int u, int v) {
        succ[u].erase(v);
        pred[v].erase(u);
    }

    std::vector<int> nodes() const {
        std::vector<int> vertices;
        for (int i = 0; i < id_counter; i++) {
            vertices.push_back(i);
        }
        return vertices;
    }

    std::vector<std::pair<int, int>> edges() const {
        std::vector<std::pair<int, int>> edges;
        for (const auto &[u, vs] : succ) {
            for (const auto &v : vs) {
                edges.push_back(std::make_pair(u, v));
            }
        }
        return edges;
    }

    /* TODO
       WARNING: does not maintain invariant
       that all edges are between 0 and N. We 
       will need to update the code to fix this.
       Probably should return a map from old to new
       vertex labels.
     */
    void delete_vertex(int u) {
        // removes u and all (v, u) edges
        vertices.erase(u);
        for (const vertex<T>& v : pred[u]) {
            succ[v].erase(u);
        }

        succ.erase(u);
        pred.erase(u);
    }

    vertex<T>& operator[](int u) {
        return vertices.at(u);
    }

    const vertex<T>& operator[](int u) const {
        return vertices.at(u);
    }

    const std::set<int>& predecessors(int u) const {
        return pred.at(u);
    }

    const std::set<int>& successors(int u) const {
        return succ.at(u);
    }

    bool contains(int u) const {
        return vertices.find(u) != vertices.end();
    }

    size_t out_degree(int u) const {
        return succ.at(u).size();
    }

    /*
     * Returns true if u is an ancestor of v in the given tree.
     */
    friend bool ancestor(const digraph<T>& tree, int u, int v) {
        if (u == v) {
            return true;
        }

        for (int w : tree.succ.at(u)) {
            if (ancestor(tree, w, v)) {
                return true;
            }
        }

        return false;
    }
};

template <class T>
std::string to_adjacency_list(const digraph<T>& G, const std::unordered_map<int, int>& vertex_map) {
    std::stringstream ss;
    for (const auto& u : G.nodes()) {
        ss << vertex_map.at(u) << " ";
        for (const auto& v : G.successors(u)) {
            ss << vertex_map.at(v) << " ";
        }
        ss << std::endl;
    }
    return ss.str();
}

/*
 * Parses an adjacency list into a directed graph object, 
 * where the vertices are read in as integers.
 */
std::pair<digraph<int>, std::unordered_map<int, int>> parse_adjacency_list(const std::string& filename) {
    digraph<int> g;

    std::ifstream file(filename);
    std::string line;
    std::unordered_map<int, int> vertex_map;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    while (std::getline(file, line)) {
        if (line.length() < 1 || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        int src, tgt;

        if (!(iss >> src)) {
            break;
        }

        if (vertex_map.find(src) == vertex_map.end()) {
            vertex_map[src] = g.add_vertex(src);
        }

        while (iss >> tgt) {
            if (vertex_map.find(tgt) == vertex_map.end()) {
                vertex_map[tgt] = g.add_vertex(tgt);
            }

            g.add_edge(vertex_map[src], vertex_map[tgt]);
        }
    }

    file.close();
    return std::make_pair(g, vertex_map);
}

/*
 * Generates a random integer in the range [a, b].
 */
int rand_int(std::ranlux48_base& gen, int a, int b) {
    std::uniform_int_distribution<int> distrib(a, b);
    return distrib(gen);
}

/*
 * Generates a random predecessor of v in the given graph, where 
 * the probability of selecting a predecessor is proportional to
 * the weight of the edge from the predecessor to v.
 */
int rand_predecessor(
        int v, const std::set<int>& predecessors, 
        const std::map<std::pair<int,int>, double> &weights, std::ranlux48_base& gen
) {
    std::vector<int> preds(predecessors.begin(), predecessors.end());
    std::vector<int> pred_weights(predecessors.size());
    for (int i = 0; i < predecessors.size(); i++) {
        pred_weights[i] = weights.at(std::make_pair(preds[i], v));
    }

    std::discrete_distribution<int> distrib(pred_weights.begin(), pred_weights.end());
    return preds[distrib(gen)];
}

/*
 * Generates a random spanning tree of the given graph using
 * Wilson's algorithm, where all edge weights are assumed to be 1.
 *
 * Returns:
 *   - Random spanning tree of G
 *   - Root of the tree
 *
 * Warning: The graph must be strongly connected.
 */
template <class T>
std::pair<digraph<T>, int> sample_random_spanning_tree(
        const digraph<T> &G, const std::map<std::pair<int,int>, double> &weights, std::ranlux48_base& gen, int root = -1
) {
    digraph<T> spanning_tree;
    for (int u : G.nodes()) {
        spanning_tree.add_vertex(G[u].data);
    }

    std::vector<int> next(G.nodes().size(), -1);
    std::vector<bool> in_tree(G.nodes().size(), false);

    if (root == -1) {
        root = rand_int(gen, 0, G.nodes().size() - 1);
    }

    in_tree[root] = true;
    for (int u : G.nodes()) {
        if (in_tree[u]) continue;

        int v = u;
        while (!in_tree[v]) {
            auto& pred = G.predecessors(v);
            if (pred.size() == 0) {
                throw std::runtime_error("Graph is not strongly connected");
            }

            next[v] = *std::next(pred.begin(), rand_predecessor(v, pred, weights, gen));
            v = next[v];
        }

        v = u;
        while (!in_tree[v]) {
            in_tree[v] = true;
            spanning_tree.add_edge(next[v], v);
            v = next[v];
        }
    }

    return std::make_pair(spanning_tree, root);
}

#endif
