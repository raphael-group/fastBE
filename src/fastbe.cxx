#include <pprint.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/ostr.h>

#include <argparse/argparse.hpp>
#include <nanobench.h>

#include <fastbe.hpp>
#include <digraph.hpp>
#include <piecewiselinearf.hpp>
#include <piecewisequadraticf.hpp>

#include <cmath>
#include <unordered_map>
#include <chrono>
#include <functional>
#include <random>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <optional>
#include <stack>
#include <tuple>

using json = nlohmann::json;

/*
 * Stores the piecewise linear functions describing 
 * $Score(id, \psi)$ for all rows of the frequency matrix
 * for a vertex in the clone tree.
 *
 * Members:
 * - id: The id of the vertex in the clone tree.
 * - valid: Whether or not the piecewise linear functions and w are up to date.
 * - w: The weights for column id of the W matrix, which is derived from F.
 * - fs: The piecewise linear functions for each row of the frequency matrix.
 */
struct clone_tree_vertex {
    int id;
    bool valid;
    std::vector<float> w;
    std::vector<PiecewiseLinearF> fs;

    clone_tree_vertex(size_t nrows, int id) : id(id), valid(false) {
        w = std::vector<float>(nrows, 0.0);
        fs = std::vector<PiecewiseLinearF>(nrows, PiecewiseLinearF({0.0}, 0.0));
    }

    clone_tree_vertex() : id(-1), valid(false) {}
};

struct clone_tree_vertex_l2 {
    int id;
    bool valid;
    std::vector<float> w;
    std::vector<PiecewiseQuadraticF> fs;

    clone_tree_vertex_l2(size_t nrows, int id) : id(id), valid(false) {
        w = std::vector<float>(nrows, 0.0);
        fs = std::vector<PiecewiseQuadraticF>(nrows, PiecewiseQuadraticF());
    }

    clone_tree_vertex_l2() : id(-1), valid(false) {}
};

/*
 * Converts a clone tree with integer vertex labels to a clone tree
 * with clone_tree_vertex vertex labels.
 *
 * To be used to avoid recomputation when calling `one_fastbe_l1` multiple 
 * times with perturbations of the same clone tree.
 */
digraph<clone_tree_vertex> convert_clone_tree(const digraph<int>& clone_tree, size_t nrows) {
    digraph<clone_tree_vertex> new_clone_tree;
    for (size_t node_id = 0; node_id < clone_tree.nodes().size(); ++node_id) {
        int clone_id = clone_tree[node_id].data;
        clone_tree_vertex vertex_data(nrows, clone_id);
        new_clone_tree.add_vertex(vertex_data);
    }

    for (auto [u, v] : clone_tree.edges()) {
        new_clone_tree.add_edge(u, v);
    }

    return new_clone_tree;
}

digraph<clone_tree_vertex_l2> convert_clone_tree_l2(const digraph<int>& clone_tree, size_t nrows) {
    digraph<clone_tree_vertex_l2> new_clone_tree;
    for (size_t node_id = 0; node_id < clone_tree.nodes().size(); ++node_id) {
        int clone_id = clone_tree[node_id].data;
        clone_tree_vertex_l2 vertex_data(nrows, clone_id);
        new_clone_tree.add_vertex(vertex_data);
    }

    for (auto [u, v] : clone_tree.edges()) {
        new_clone_tree.add_edge(u, v);
    }

    return new_clone_tree;
}

/*
 * Invalidates all vertices on the path from u to root.
 */ 
void invalidate(digraph<clone_tree_vertex>& tree, int u, int root) {
    tree[u].data.valid = false;

    if (u == root) return;

    int parent = *tree.predecessors(u).begin();
    invalidate(tree, parent, root);
}

/*
 * Given a frequency matrix $F$ and a clone tree $T$, this function
 * finds the minimizing value of $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_2^2$$ 
 * over all usage matrices in $\mathcal{O}(mn^2)$ using dynamic programming.
 *
 * This function avoids recomputation, by only updating the 
 * piecewise linear functions of the clone tree at vertices 
 * whose subtrees have been modified, as specified by the 
 * `valid` flag in the `clone_tree_vertex` struct.
 *
 * Input: 
 *  - clone_tree: A clone tree represented as a digraph.
 *  - vertex_map: A map from the columns of the frequency matrix to the 
 *    vertices of the clone tree.
 *  - F: A frequency matrix represented as a 2D vector.
*/
float one_fastbe_l2(
    digraph<clone_tree_vertex_l2>& clone_tree, 
    const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F, 
    int root
) {
    size_t nrows = F.size();

    std::stack<int> call_stack; // contains vertices to visit in column coordinates
    call_stack.push(root);
    while(!call_stack.empty()) {
        int i = call_stack.top();
        call_stack.pop();

        if (clone_tree[vertex_map.at(i)].data.valid) continue;

        // If leaf, set piecewise linear function and return.
        if (clone_tree.out_degree(vertex_map.at(i)) == 0) {
            for (size_t j = 0; j < nrows; ++j) {
                clone_tree[vertex_map.at(i)].data.fs[j] = PiecewiseQuadraticF(F[j][i]);
            }

            clone_tree[vertex_map.at(i)].data.valid = true;
            continue;
        }

        // Recurse at children. 
        bool all_children_valid = true; 
        for (auto k : clone_tree.successors(vertex_map.at(i))) {
            if (!clone_tree[k].data.valid) {
                if (all_children_valid) {
                    call_stack.push(i);
                }

                call_stack.push(clone_tree[k].data.id);
                all_children_valid = false;
            }
        }

        if (!all_children_valid) continue;

        for (size_t j = 0; j < nrows; ++j) {
            PiecewiseQuadraticF g;
            for (auto k : clone_tree.successors(vertex_map.at(i))) {
                PiecewiseQuadraticF f = clone_tree[k].data.fs[j]; // copy!
                g = g + f; 
            }

            clone_tree[vertex_map.at(i)].data.fs[j] = g.update_representation(F[j][i]);
        }

        clone_tree[vertex_map.at(i)].data.valid = true;
    }

    float obj = 0;
    for (size_t j = 0; j < nrows; ++j) {
        PiecewiseQuadraticF f = clone_tree[vertex_map.at(root)].data.fs[j];
        std::vector<double> breakpoints = f.breakpoints;
        std::vector<double> values = f.evaluate_at_breakpoints();
        double obj_inc = f.f0;
        for (size_t i = 0; i < breakpoints.size(); ++i) {
            if (breakpoints[i] < 0.0) continue;
            obj_inc = std::max(obj_inc, values[i] - breakpoints[i]);
        }

        obj += obj_inc;
    }

    return obj;
}

/*
 * Given a frequency matrix $F$ and a clone tree $T$, this function
 * finds the minimizing value of $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$ 
 * over all usage matrices in $\mathcal{O}(mnd)$ using dynamic programming.
 *
 * This function avoids recomputation, by only updating the 
 * piecewise linear functions of the clone tree at vertices 
 * whose subtrees have been modified, as specified by the 
 * `valid` flag in the `clone_tree_vertex` struct.
 *
 * Input: 
 *  - clone_tree: A clone tree represented as a digraph.
 *  - vertex_map: A map from the columns of the frequency matrix to the 
 *    vertices of the clone tree.
 *  - F: A frequency matrix represented as a 2D vector.
*/
float one_fastbe_l1(
    digraph<clone_tree_vertex>& clone_tree, 
    const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F, 
    int root
) {
    size_t nrows = F.size();

    std::stack<int> call_stack; // contains vertices to visit in column coordinates
    call_stack.push(root);
    while(!call_stack.empty()) {
        int i = call_stack.top();
        call_stack.pop();

        if (clone_tree[vertex_map.at(i)].data.valid) continue;

        std::vector<float>& w_ref = clone_tree[vertex_map.at(i)].data.w;

        /* if not valid, first set w vector. */
        for (size_t j = 0; j < nrows; ++j) {
            w_ref[j] = F[j][i];
            for (auto k : clone_tree.successors(vertex_map.at(i))) {
                w_ref[j] -= F[j][clone_tree[k].data.id];
            }
        }

        /* If leaf, set piecewise linear function and return. */
        if (clone_tree.out_degree(vertex_map.at(i)) == 0) {
            for (size_t j = 0; j < nrows; ++j) {
                clone_tree[vertex_map.at(i)].data.fs[j] = PiecewiseLinearF({w_ref[j]}, 0);
            }

            clone_tree[vertex_map.at(i)].data.valid = true;
            continue;
        }

        /* Recurse at children. */
        bool all_children_valid = true; 
        for (auto k : clone_tree.successors(vertex_map.at(i))) {
            if (!clone_tree[k].data.valid) {
                if (all_children_valid) {
                    call_stack.push(i);
                }

                call_stack.push(clone_tree[k].data.id);
                all_children_valid = false;
            }
        }

        if (!all_children_valid) continue;

        for (size_t j = 0; j < nrows; ++j) {
            PiecewiseLinearF g_out({w_ref[j]}, 0);
            for (auto k : clone_tree.successors(vertex_map.at(i))) {
                PiecewiseLinearF f = clone_tree[k].data.fs[j]; // copy!
                f.compute_minimizer();
                g_out.addInPlace(std::move(f));
            }

            clone_tree[vertex_map.at(i)].data.fs[j] = g_out;
        }

        clone_tree[vertex_map.at(i)].data.valid = true;
    }

    float obj = 0;
    for (size_t j = 0; j < nrows; ++j) {
        PiecewiseLinearF f = clone_tree[vertex_map.at(root)].data.fs[j];
        f.compute_minimizer();
        f.addInPlace(PiecewiseLinearF({1 - F[j][root]}, 0));
        obj += f.minimizer();
    }

    return -1 * obj;
}

/* 
 * Computes the score of all 2^{|child(v)|} ways to place vertex u 
 * as a child of vertex v when allowing the children of v to be 
 * placed as children of u.
 */
std::vector<std::pair<float, long long>> score_placements(
    digraph<clone_tree_vertex> partial_tree,
    const std::unordered_map<int, int> &vertex_map,
    const std::vector<std::vector<float>>& F, 
    int root_index, int root, int u, int v
) {
    std::vector<int> v_children = partial_tree.successors(v);
    std::vector<int> children_parents(v_children.size(), v);

    partial_tree.add_edge(v, u);

    std::vector<std::pair<float, long long>> placements;
    for (long long i=0; i<(1<<v_children.size()); i++) {
        for (size_t j = 0; j < v_children.size(); j++) {
            partial_tree.remove_edge(children_parents[j], v_children[j]);
        }

        for (size_t j=0; j<v_children.size(); j++) {
            if (i & (1<<j)) {
                partial_tree.add_edge(v, v_children[j]);
                children_parents[j] = v;
            } else {
                partial_tree.add_edge(u, v_children[j]);
                children_parents[j] = u;
            }
        }

        invalidate(partial_tree, u, root_index);
        float score = one_fastbe_l1(partial_tree, vertex_map, F, root);
        placements.push_back({score, i});
    }

    return placements;
}

/* 
 * Performs (forward) beam search to infer a clone tree from a frequency matrix.
 */
std::pair<digraph<clone_tree_vertex>, std::unordered_map<int, int>> forward_beam_search(
    const std::vector<std::vector<float>>& F, 
    int root,
    std::vector<int> clone_order,
    size_t beam_width,
    unsigned int num_threads
) {
    std::unordered_map<int, int> vertex_map; // maps from the clone id (column idx) to the vertex id in the clone tree
    std::vector<digraph<clone_tree_vertex>> partial_trees;

    { 
        digraph<clone_tree_vertex> initial_tree; // starts with all clones and no edges
        for (auto clone : clone_order) {
            vertex_map[clone] = initial_tree.add_vertex(clone_tree_vertex(F.size(), clone));
        }
        partial_trees.push_back(initial_tree);
    }

    int root_index = vertex_map[root];

    int counter = 1;
    for (size_t j = 0; j < clone_order.size(); ++j) {
        auto clone = clone_order[j];

        spdlog::info("Adding clone {} ({}/{}) to the partially constructed trees", clone, counter, clone_order.size());
        counter++;

        if (clone == root) continue;

        std::mutex proposed_trees_mutex;
        std::vector<std::tuple<float, size_t, int, long long>> proposed_trees; // (score, partial tree idx, parent vertex, child placement)
        int u = vertex_map[clone];

        for (size_t i = 0; i < partial_trees.size(); ++i) {
            std::vector<std::thread> threads;
            for (size_t k = 0; k < num_threads; ++k) {
                threads.push_back(std::thread([&](size_t i, size_t k) {
                    for (size_t l = k; l < j; l += num_threads) {
                        int v = vertex_map[clone_order[l]];
                        auto placements = score_placements(partial_trees[i], vertex_map, F, root_index, root, u, v);
                        for (auto [score, placement] : placements) {
                            std::lock_guard<std::mutex> lock(proposed_trees_mutex);
                            proposed_trees.push_back(std::make_tuple(score, i, v, placement));
                        }
                    }
                }, i, k));
            }

            for (auto& thread : threads) {
                thread.join();
            }
        }

        std::sort(proposed_trees.begin(), proposed_trees.end(), [&](auto& a, auto& b) {
            return std::get<0>(a) < std::get<0>(b);
        });
        
        std::vector<digraph<clone_tree_vertex>> new_partial_trees;
        for (auto [score, i, v, placement] : proposed_trees) {
            auto tree = partial_trees[i];
            if (new_partial_trees.size() == beam_width) break;
            
            std::vector<int> v_children = tree.successors(v);
           
            // remove all edges from v
            for (auto w : v_children) {
                tree.remove_edge(v, w);
            }

            // add edges according to placement
            tree.add_edge(v, u);
            for (size_t j = 0; j < v_children.size(); j++) {
                if (placement & (1<<j)) {
                    tree.add_edge(v, v_children[j]);
                } else {
                    tree.add_edge(u, v_children[j]);
                }
            }
    
            invalidate(tree, u, root_index);
            new_partial_trees.push_back(tree);
        }

        partial_trees = new_partial_trees;
    }

    std::sort(partial_trees.begin(), partial_trees.end(), [&](auto& a, auto& b) {
        return one_fastbe_l1(a, vertex_map, F, root) < one_fastbe_l1(b, vertex_map, F, root);
    });

    spdlog::info("Best tree objective is {}", one_fastbe_l1(partial_trees[0], vertex_map, F, root));
    spdlog::info("Worst tree objective is {}", one_fastbe_l1(partial_trees[partial_trees.size() - 1], vertex_map, F, root));
    
    return {partial_trees[0], vertex_map};
}

/* 
 * This function parses a frequency matrix from a text file and returns it as a 2D vector.
 *
 * Input: 
 *  - filename: The name of the file that contains the frequency matrix.
 *              Each line in the file should correspond to a row in the matrix.
 *              Within a line, numbers should be space-separated and represent the values 
 *              in the columns of that row.
 * 
 * Returns a 2D vector (vector of vectors) that represents the frequency matrix. Each inner 
 * vector represents a row of the matrix. Throws a std::runtime_error if the file can't be opened
 * or if the file does not represent a valid matrix.
*/
std::vector<std::vector<float>> parse_frequency_matrix(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the file.");
    }

    std::vector<std::vector<float>> matrix;
    std::string line;
    size_t num_cols = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> row;
        float value;

        while (iss >> value) {
            row.push_back(value);
        }

        if (matrix.empty()) {
            num_cols = row.size();
        } else if (row.size() != num_cols) {
            throw std::runtime_error("The file does not represent a matrix. Number of columns is not consistent across rows.");
        }

        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

    
enum class loss_type {L1, L2};

float sum(const std::vector<std::vector<float>>& frequency_matrix, const std::vector<int>& columns, int i) {
    float sum = 0;
    for (auto col : columns) {
        sum += frequency_matrix[i][col];
    }
    return sum;
}

float median(const std::vector<std::vector<float>>& frequency_matrix, const std::vector<int>& columns, int i) {
    std::vector<float> values;
    for (auto col : columns) {
        values.push_back(frequency_matrix[i][col]);
    }

    std::sort(values.begin(), values.end());
    size_t n = values.size();
    if (n % 2 == 0) {
        return (values[n/2 - 1] + values[n/2]) / 2;
    } else {
        return values[n/2];
    }
}

/*
 * Computes the score of a cluster of clones (columns) in the 
 * frequency matrix.
 */
float compute_component_score(
    loss_type loss,
    const std::vector<std::vector<float>>& frequency_matrix,
    const std::vector<int>& columns
) {
    int m = frequency_matrix.size();

    // compute optimal center
    std::vector<float> center(m, 0.0);
    for (int i = 0; i < m; ++i) {
        if (loss == loss_type::L2) {
            center[i] = sum(frequency_matrix, columns, i) / columns.size();
        } else {
            center[i] = median(frequency_matrix, columns, i);
        }
    }

    double obj = 0;
    for (int i = 0; i < m; i++) {
        for (auto j : columns) {
            if (loss == loss_type::L2) {
                obj += (frequency_matrix[i][j] - center[i]) * (frequency_matrix[i][j] - center[i]);
            } else {
                obj += std::abs(frequency_matrix[i][j] - center[i]);
            }
        }
    }

    return obj;
}

/* 
 * Returns the connected components of an induced subgraph 
 * in the clone tree when viewed as an undirected graph.
 */
std::vector<std::vector<int>> connected_components(
    const digraph<int>& clone_tree,
    const std::vector<int>& subgraph
) {
    std::set<int> unvisited(subgraph.begin(), subgraph.end());
    std::vector<std::vector<int>> components;

    while (!unvisited.empty()) {
        std::queue<int> q;
        q.push(*unvisited.begin());
        std::vector<int> component;
        while (!q.empty()) {
            int node = q.front();
            q.pop();
            component.push_back(node);
            unvisited.erase(node);

            for (auto neighbor : clone_tree.successors(node)) {
                if (unvisited.find(neighbor) != unvisited.end()) { // if not visited
                    q.push(neighbor);
                }
            }

            for (auto neighbor : clone_tree.predecessors(node)) {
                if (unvisited.find(neighbor) != unvisited.end()) {
                    q.push(neighbor);
                }
            }
        }
        components.push_back(component);
    }

    return components;
}

/*
 * Given a clone tree $T$ and a frequency matrix $F$, this function,
 * and a parameter $k$, this function computes clustering of the columns 
 * of the frequency matrix into $1, \ldots, k$ clusters which is
 * consistent with the clone tree $T$.
 *
 * The function is not guaranteed to return exactly $k$ clusters if
 * optimal clustering is possible in with fewer clusters.
 */
std::pair<std::vector<std::pair<float, std::vector<int>>>, std::vector<float>> divisive_clustering(
    loss_type loss,
    const std::vector<std::vector<float>>& frequency_matrix,
    digraph<int>& clone_tree,
    size_t k
) {
    if (k == 0) return {{}, {}};

    std::set<std::pair<float, std::vector<int>>> clustering;

    // create initial component consisting of all columns 
    std::vector<int> initial_component = clone_tree.nodes();
    std::vector<int> initial_columns(initial_component.size());
    for (size_t i = 0; i < initial_component.size(); ++i) {
        initial_columns[i] = clone_tree[initial_component[i]].data;
    }

    float initial_obj = compute_component_score(loss, frequency_matrix, initial_columns);
    clustering.insert({-initial_obj, initial_columns});

    // in theory, can implement in O(mnk) time, this implementation is O(mn^2klog(n))
    std::vector<float> obj_values(k, 0.000);
    obj_values[0] = initial_obj;
    while(clustering.size() < k) {
        spdlog::info("Objective value with {} clusters: {}", clustering.size(), obj_values[clustering.size() - 1]);

        auto [obj, component] = *clustering.begin();
        if (component.size() == 1) break;

        // TODO: a little slow, make component a set
        std::vector<std::pair<int, int>> component_edges;
        for (auto edge : clone_tree.edges()) {
            if (std::find(component.begin(), component.end(), edge.first) != component.end() &&
                std::find(component.begin(), component.end(), edge.second) != component.end()) {
                component_edges.push_back(edge);
            }
        }

        // TODO: can improve by not recomputing objective
        // for each split. instead, walk through tree maintaining
        // the current objective and the current split and update 
        // it as we go in O(1) time using median/mean update formulas.
        //
        // this will allow us to find the best split in O(n) time. though
        // arguably, this is simpler.
        float min_obj1 = std::numeric_limits<float>::max();
        float min_obj2 = std::numeric_limits<float>::max();
        float min_obj  = std::numeric_limits<float>::max();
        std::vector<int> sub_component1, sub_component2;
        for (auto [u, v] : component_edges) {
            clone_tree.remove_edge(u, v);
            auto ccs = connected_components(clone_tree, component);
            clone_tree.add_edge(u, v);

            assert(ccs.size() == 2);

            auto c1 = ccs[0];
            auto c2 = ccs[1];

            std::vector<int> c1_ids(c1.size());
            std::vector<int> c2_ids(c2.size());
            for (size_t i = 0; i < c1.size(); ++i) {
                c1_ids[i] = clone_tree[c1[i]].data;
            }
            for (size_t i = 0; i < c2.size(); ++i) {
                c2_ids[i] = clone_tree[c2[i]].data;
            }

            float obj1 = compute_component_score(loss, frequency_matrix, c1_ids);
            float obj2 = compute_component_score(loss, frequency_matrix, c2_ids);

            if (obj1 + obj2 < min_obj) {
                min_obj = obj1 + obj2;
                min_obj1 = obj1;
                min_obj2 = obj2;
                sub_component1 = c1;
                sub_component2 = c2;
            }
        }

        clustering.erase(clustering.begin());
        clustering.insert({-min_obj1, sub_component1});
        clustering.insert({-min_obj2, sub_component2});

        // update list of objective values
        float current_obj = 0;
        for (auto [obj, component] : clustering) {
            current_obj += obj;
        }

        obj_values[clustering.size() - 1] = -current_obj;
    }

    spdlog::info("Objective value with {} clusters: {}", clustering.size(), obj_values[clustering.size() - 1]);

    std::vector<std::pair<float, std::vector<int>>> new_clustering;
    for (auto [obj, component] : clustering) {
        new_clustering.push_back({-obj, component});
    }

    return {new_clustering, obj_values};
}

void perform_cluster(const argparse::ArgumentParser &cluster) {
    size_t num_clusters = cluster.get<size_t>("num_clusters");

    auto [clone_tree_int, vertex_map] = parse_adjacency_list(cluster.get<std::string>("clone_tree"));
    std::vector<std::vector<float>> frequency_matrix = parse_frequency_matrix(cluster.get<std::string>("frequency_matrix"));

    loss_type loss = loss_type::L1;
    if (cluster.get<std::string>("loss") == "L2") {
        loss = loss_type::L2;
    } else if (cluster.get<std::string>("loss") != "L1") {
        throw std::runtime_error("Invalid loss function specified.");
    }

    auto start = std::chrono::high_resolution_clock::now();
    spdlog::info("Performing divisive clustering...");
    auto [clustering, obj_values] = divisive_clustering(
        loss,
        frequency_matrix, 
        clone_tree_int, 
        num_clusters
    );

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    std::ofstream adj_output(cluster.get<std::string>("output") + "_clustering.csv");
    adj_output << "mutation,clone" << std::endl;
    for (size_t i = 0; i < clustering.size(); ++i) {
        for (auto vert : clustering[i].second) {
            adj_output << clone_tree_int[vert].data << "," << i << std::endl;
        }
    }

    json res;
    res["time (ns)"] = duration.count();
    res["time (s)"]  = duration.count() / 1e9;
    res["num_clusters"] = clustering.size();
    res["objective_values"] = obj_values;
    res["loss_function"] = loss == loss_type::L1 ? "L1" : "L2";

    std::ofstream json_output_stats(cluster.get<std::string>("output") + "_clustering_results.json");
    json_output_stats << res.dump(4) << std::endl;
}

void perform_regression(const argparse::ArgumentParser &regress) {
    auto [clone_tree_int, vertex_map] = parse_adjacency_list(regress.get<std::string>("clone_tree"));
    std::vector<std::vector<float>> frequency_matrix = parse_frequency_matrix(regress.get<std::string>("frequency_matrix"));
    auto clone_tree = convert_clone_tree_l2(clone_tree_int, frequency_matrix.size());

    int root = regress.get<int>("assigned_root");
    auto start = std::chrono::high_resolution_clock::now();
    float obj = one_fastbe_l2(clone_tree, vertex_map, frequency_matrix, root);

    for (size_t i = 0; i < regress.get<size_t>("num_reps") - 1; ++i) {
        for (auto vertex : clone_tree.nodes()) {
            clone_tree[vertex].data.valid = false;
        }

        one_fastbe_l2(clone_tree, vertex_map, frequency_matrix, root);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    json output;
    output["time (ns)"] = duration.count();
    output["time (s)"]  = duration.count() / 1e9;
    output["objective_value"] = obj;

    std::ofstream output_file(regress.get<std::string>("output") + "_results.json");
    output_file << output.dump(4) << std::endl;
    output_file.close();

    return;
}

void perform_search(const argparse::ArgumentParser &search) {
    std::vector<std::vector<float>> frequency_matrix = parse_frequency_matrix(search.get<std::string>("frequency_matrix"));

    size_t nrows = frequency_matrix.size();
    size_t ncols = frequency_matrix[0].size();

    unsigned int num_threads = search.get<unsigned int>("threads");

    std::vector<int> clone_order(ncols); // map from order to clone
    for (size_t i = 0; i < ncols; ++i) {
        clone_order[i] = i;
    }

    // use F-sum ordering to order clones
    std::sort(clone_order.begin(), clone_order.end(), [&](int i, int j) {
        float sum_i = 0.0;
        float sum_j = 0.0;
        for (size_t k = 0; k < nrows; ++k) {
            sum_i += frequency_matrix[k][i];
            sum_j += frequency_matrix[k][j];
        }
        return sum_i > sum_j;
    });

    int root = search.get<int>("assigned_root");
    auto root_it = std::find(clone_order.begin(), clone_order.end(), root);
    clone_order.erase(root_it);
    clone_order.insert(clone_order.begin(), root);
    
    int beam_width = search.get<int>("beam_width");
    if (beam_width == -1) {
        if (ncols <= 10) {
            beam_width = 1000;
        } else if (ncols <= 25) {
            beam_width = 500;
        } else if (ncols <= 50) {
            beam_width = 250;
        } else if (ncols <= 100) {
            beam_width = 100;
        } else if (ncols <= 200) {
            beam_width = 25;
        } else if (ncols <= 500) {
            beam_width = 5;
        } else {
            beam_width = 1;
        }

        spdlog::info("Beam width not specified. Using default beam width of {} for {} clones.", beam_width, ncols);
    }

    spdlog::info("Performing beam search to find tree(s)...");
    auto [clone_tree, vmap] = forward_beam_search(
        frequency_matrix, 
        root, 
        clone_order, 
        beam_width,
        num_threads
    );

    float obj = one_fastbe_l1(clone_tree, vmap, frequency_matrix, root);
    spdlog::info("Objective value: {}", obj);

    std::string adjacency_list = to_adjacency_list(clone_tree, vmap);

    std::ofstream adj_output(search.get<std::string>("output") + "_tree.txt");
    adj_output << adjacency_list;

    json output_stats;
    output_stats["objective_value"] = obj;
    std::ofstream json_output_stats(search.get<std::string>("output") + "_results.json");
    json_output_stats << output_stats.dump(4) << std::endl;

    return;
}

int main(int argc, char *argv[]) {
    auto console_logger = spdlog::stdout_color_mt("fastbe");
    spdlog::set_default_logger(console_logger);

    auto error_logger = spdlog::stderr_color_mt("error");

    argparse::ArgumentParser program(
        "fastbe",
        std::to_string(FASTBE_VERSION_MAJOR) + "." + std::to_string(FASTBE_VERSION_MINOR)
    );

    argparse::ArgumentParser regress(
        "regress"
    );

    argparse::ArgumentParser search(
        "search"
    );

    argparse::ArgumentParser cluster(
        "cluster"
    );


    regress.add_description("regresses a clone tree onto a frequency matrix.");
    search.add_description("searches for a clone tree that best fits a frequency matrix.");
    cluster.add_description("clusters the mutations (columns) of a frequency matrix.");

    regress.add_argument("clone_tree")
           .help("adjacency list of the clone tree");

    regress.add_argument("frequency_matrix")
           .help("TXT file containing the frequency matrix");

    regress.add_argument("-o", "--output")
           .help("prefix of the output files")
           .required();

    regress.add_argument("-n", "--num_reps")
           .help("number of times to repeat the regression for benchmarking")
           .default_value((size_t) 1)
           .scan<'u', size_t>();

    regress.add_argument("-f", "--assigned_root")
          .help("assigned root")
          .default_value(0)
          .scan<'d', int>();

    search.add_argument("frequency_matrix")
          .help("TXT file containing the frequency matrix");

    search.add_argument("-o", "--output")
          .help("prefix of the output files")
          .required();

    search.add_argument("-t", "--threads")
          .help("number of threads to use")
          .default_value(std::thread::hardware_concurrency())
          .scan<'u', unsigned int>();

    search.add_argument("-f", "--assigned_root")
          .help("assigned root")
          .default_value(0)
          .scan<'d', int>();

    search.add_argument("-b", "--beam_width")
          .help("beam width")
          .default_value(-1)
          .scan<'d', int>();

    cluster.add_argument("clone_tree")
           .help("adjacency list of the clone tree");

    cluster.add_argument("frequency_matrix")
           .help("TXT file containing the frequency matrix");

    cluster.add_argument("-k", "--num_clusters")
           .help("number of clusters")
           .required()
           .scan<'u', size_t>();

    cluster.add_argument("-o", "--output")
           .help("prefix of the output files")
           .required();

    cluster.add_argument("-l", "--loss")
           .help("loss function to use for clustering")
           .default_value("L1")
           .choices("L1", "L2");

    program.add_subparser(search);
    program.add_subparser(regress);
    program.add_subparser(cluster);
    
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;

        if (program.is_subcommand_used(regress)) {
            std::cerr << regress;
        } else if (program.is_subcommand_used(search)) {
            std::cerr << search;
        } else if (program.is_subcommand_used(cluster)) {
            std::cerr << cluster;
        } else {
            std::cerr << program;
        }

        std::exit(1);
    }

    if (program.is_subcommand_used(regress)) {
        perform_regression(regress);
    } else if (program.is_subcommand_used(search)) {
        perform_search(search);
    } else if (program.is_subcommand_used(cluster)) {
        perform_cluster(cluster);
    } else {
        std::cerr << program;
    }

    return 0;
}
