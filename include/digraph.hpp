#ifndef _DIGRAPH_H
#define _DIGRAPH_H

#include <utility>
#include <vector>
#include <set>
#include <map>
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

    std::map<int, std::set<int>> succ;
    std::map<int, std::set<int>> pred;

    std::map<int, vertex<T>> vertices;
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
};

/*
 * Parses an adjacency list into a directed graph object, 
 * where the vertices are read in as integers.
 */
std::pair<digraph<int>, std::map<int, int>> parse_adjacency_list(const std::string& filename) {
    digraph<int> g;

    std::ifstream file(filename);
    std::string line;
    std::map<int, int> vertex_map;

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

#endif
