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

#include <cmath>
#include <unordered_map>
#include <random>
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

/*
 * Converts a clone tree with integer vertex labels to a clone tree
 * with clone_tree_vertex vertex labels.
 *
 * To be used to avoid recomputation when calling `one_fastbe` multiple 
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

/*
 * Prunes the subtree rooted at vertex u and regrafts as a child of 
 * vertex v. This function assumes that u is not the root of the tree.
 */
void subtree_prune_and_regraft(digraph<clone_tree_vertex>& tree, int u, int v, int root) {
    int parent = *tree.predecessors(u).begin();
    tree.remove_edge(parent, u);
    tree.add_edge(v, u);

    std::stack<int> s;
    s.push(u);
    s.push(parent);

    // walk from vertices u and v to root and update valid flags
    while (!s.empty()) {
        int node_id = s.top();
        s.pop();

        tree[node_id].data.valid = false;

        if (node_id == root) continue;

        int parent = *tree.predecessors(node_id).begin();
        s.push(parent);
    }
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
float one_fastbe(
    digraph<clone_tree_vertex>& clone_tree, 
    const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F, 
    int root
) {
    size_t nrows = F.size();

    std::stack<int> call_stack;
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
        PiecewiseLinearF f = clone_tree[root].data.fs[j];
        f.compute_minimizer();
        f.addInPlace(PiecewiseLinearF({1 - F[j][root]}, 0));
        obj += f.minimizer();
    }

    return -1 * obj;
}

/*
 * Given a frequency matrix $F$ and a clone tree $T$, this function
 * finds the minimizing value of $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$ 
 * over all usage matrices in $\mathcal{O}(mn^2)$ using dynamic programming.
 *
 * Input: 
 *  - clone_tree: A clone tree represented as a digraph.
 *  - vertex_map: A map from the columns of the frequency matrix to the 
 *    vertices of the clone tree.
 *  - columns: A vector of column indices to consider.
 *  - F: A frequency matrix represented as a 2D vector.
 *
 * Requirement:
 *  - The selected `columns` must coincide with the vertices of the clone tree.
*/
float one_fastbe(
    const digraph<int>& clone_tree, const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F, const std::vector<int>& columns,
    int root
) {
    size_t nrows = F.size();

    std::vector<std::vector<float>> W(nrows, std::vector<float>(columns.size(), 0.0));

    for (size_t j = 0; j < nrows; ++j) {
        for (size_t i = 0; i < columns.size(); ++i) {
            W[j][i] = F[j][columns[i]];
            for (auto k : clone_tree.successors(vertex_map.at(columns[i]))) {
                W[j][i] -= F[j][clone_tree[k].data];
            }
        }
    }

    std::unordered_map<int, int> columns_inv(columns.size());
    for (size_t i = 0; i < columns.size(); ++i) {
        columns_inv[columns[i]] = i;
    }

    // i corresponds to column we are at in F!
    std::function<PiecewiseLinearF(size_t, size_t)> one_fastbe_recursive = [&](size_t j, size_t i) {
        if (clone_tree.out_degree(vertex_map.at(i)) == 0) {
            return PiecewiseLinearF({W[j][columns_inv[i]]}, 0); 
        }

        PiecewiseLinearF g_out({W[j][columns_inv[i]]}, 0);
        for (auto k : clone_tree.successors(vertex_map.at(i))) {
            auto f = one_fastbe_recursive(j, clone_tree[k].data);
            f.compute_minimizer();
            g_out.addInPlace(std::move(f));
        }

        return g_out;
    };

    float obj = 0;
    for (size_t j = 0; j < nrows; ++j) {
        auto f = one_fastbe_recursive(j, root);
        f.compute_minimizer();
        f.addInPlace(PiecewiseLinearF({1 - F[j][0]}, 0));
        obj += f.minimizer();
    }

    return -1 * obj;
}

/*
 * Computes the matrix for the total violation of the sum condition 
 * to the frequency matrix $F$ and the clone tree $T$. In particular,
 *    $$A[i, j] = \sum_{k\in C(j)} F_{i,k}$$
 * which allows us to compute (and update) the total violation of the
 * sum condition in $\mathcal{O}(nm)$ time.
 */
template <typename T>
std::vector<std::vector<float>> compute_sum_violation_matrix(
    const digraph<T>& clone_tree, const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F_transpose
) {
    size_t nrows = F_transpose.size();
    size_t ncols = F_transpose[0].size();

    std::vector<std::vector<float>> A(nrows, std::vector<float>(ncols, 0.0));

    for (size_t j = 0; j < nrows; ++j) {
        for (size_t i = 0; i < ncols; ++i) {
            float A_ji = 0;
            for (auto k : clone_tree.successors(vertex_map.at(j))) {
                A_ji += F_transpose[clone_tree[k].data.id][i];
            }
            A[j][i] = A_ji;
        }
    }

    return A;
}

/* 
 * Performs hill climbing search to find the clone tree $T$ minimizing
 * the total violation of the sum condition with respect to the frequency
 * matrix $F$ by using subtree prune and regraft operations. 
 *
 * Stops when no SPR operation reduces the error.
 *
 * Input:
 * - clone_tree: A clone tree represented as a digraph.
 * - vertex_map: A map from the columns of the frequency matrix to the
 *   vertices of the clone tree.
 * - F: A frequency matrix represented as a 2D vector.
 * - root: The root of the clone tree.
 */
void deterministic_total_violation_hill_climb(
    digraph<clone_tree_vertex>& clone_tree, const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F, int root, int max_iterations, size_t sample_id, 
    int progress_interval = 10
) {
    size_t ncols = F[0].size(), nrows = F.size();

    if (max_iterations == -1) {
        max_iterations = std::numeric_limits<int>::max();
    }

    std::vector<std::vector<float>> F_transpose(ncols, std::vector<float>(nrows, 0.0));
    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            F_transpose[j][i] = F[i][j];
        }
    }

    auto A = compute_sum_violation_matrix(clone_tree, vertex_map, F_transpose);
    float total_violation = 0;
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t j = 0; j < nrows; ++j) {
            total_violation += std::max(A[i][j] - F_transpose[i][j], 0.0f);
        }
    }

    /* 
     * Invariant(s): 
     *    - total_violation: at the beginning and end of every iteration is equal to the
     *    total violation of the sum condition.
     */
    int count = 0;
    while (true) {
        if (count % progress_interval == 0) {
            spdlog::info("Sample ID {}: iteration {}, total violation: {}", sample_id, count, total_violation);
        }

        /*
         * run DFS to compute a pre-order traversal of the tree
         * to enable constant time ancestor queries
         */
        std::stack<int> call_stack;
        call_stack.push(root);
        std::vector<int> order(clone_tree.nodes().size());
        int order_index = 0;
        while(!call_stack.empty()) {
            int i = call_stack.top();
            call_stack.pop();

            order[i] = order_index;
            order_index++;

            for (auto k : clone_tree.successors(i)) {
                call_stack.push(k);
            }
        }

        /* 
         * Go over all possible moves, greedily choosing the move that maximally 
         * reduces the total violation score. Maintains the invariant that
         * `total_violation` is equal to the total violation of the sum condition
         * for the current tree, even after the move is made.
         */
        float best_total_violation = total_violation;
        std::pair<int, int> best_move = {-1, -1};
        for (auto u : clone_tree.nodes()) {
            for (auto v : clone_tree.nodes()) {
                // u is an ancestor of v if and only if order[u] <= order[v]
                if (u == v || u == root || order[u] <= order[v]) {
                    continue; 
                }

                int p = *clone_tree.predecessors(u).begin();
                int p_column = clone_tree[p].data.id;
                int u_column = clone_tree[u].data.id;
                int v_column = clone_tree[v].data.id;

                float diff = 0;
                for (size_t j = 0; j < A[0].size(); ++j) {
                    // simulating the SPR move by putting u as a child of v
                    // i.e. A[j][p] <- A[j][p] - F[j][u]
                    float A_jp = A[p_column][j] - F_transpose[u_column][j];
                    diff += std::max(A_jp - F_transpose[p_column][j], 0.0f) - std::max(A[p_column][j] - F_transpose[p_column][j], 0.0f);
                }

                for (size_t j = 0; j < A[0].size(); ++j) {
                    // simulating the SPR move by putting u as a child of v
                    // i.e. A[j][v] <- A[j][v] + F[j][u]
                    float A_jv = A[v_column][j] + F_transpose[u_column][j];
                    diff += std::max(A_jv - F_transpose[v_column][j], 0.0f) - std::max(A[v_column][j] - F_transpose[v_column][j], 0.0f);
                }

                if (diff + total_violation < best_total_violation - 1e-4) {
                    best_total_violation = diff + total_violation;
                    best_move = {u, v};
                }
            }
        }

        if (best_move.first == -1) {
            spdlog::info("Sample ID {} found local minimum -- no SPR moves can improve the score.", sample_id);
            break;
        }

        auto [u, v] = best_move;
        subtree_prune_and_regraft(clone_tree, u, v, root);

        // update the total violation matrix, asymptotically acceptable
        // since we already spend O(n^2m) time going over all possible moves 
        A = compute_sum_violation_matrix(clone_tree, vertex_map, F_transpose);
        total_violation = 0;
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[0].size(); ++j) {
                total_violation += std::max(A[i][j] - F_transpose[i][j], 0.0f);
            }
        }

        count++;
    }

    spdlog::info("Sample ID {} hill climbing completed @ {} iterations", sample_id, count);
}

/* 
 * Performs hill climbing search to find clone tree $T$ and usage 
 * matrix $U$ minimizing  the reconstruction error $||F - UB_{T}||_1$ 
 * of the clone tree by using subtree prune and  regraft  operations.  
 *
 * Stops when no SPR operation improves the reconstruction error.
 *
 * Input:
 * - clone_tree: A clone tree represented as a digraph.
 * - vertex_map: A map from the columns of the frequency matrix to the
 *   vertices of the clone tree.
 * - F: A frequency matrix represented as a 2D vector.
 * - root: The root of the clone tree.
 */
void deterministic_hill_climb(
    digraph<clone_tree_vertex>& clone_tree, const std::unordered_map<int, int>& vertex_map, 
    const std::vector<std::vector<float>>& F, int root, int max_iterations, size_t sample_id, 
    int progress_interval = 10
) {
    size_t ncols = F[0].size(), nrows = F.size();

    if (max_iterations == -1) {
        max_iterations = std::numeric_limits<int>::max();
    }

    std::vector<std::vector<float>> F_transpose(ncols, std::vector<float>(nrows, 0.0));
    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            F_transpose[j][i] = F[i][j];
        }
    }

    auto A = compute_sum_violation_matrix(clone_tree, vertex_map, F_transpose);
    float total_violation = 0;
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t j = 0; j < nrows; ++j) {
            total_violation += std::max(A[i][j] - F_transpose[i][j], 0.0f);
        }
    }

    /* 
     * Invariant(s): 
     *    - total_violation: at the beginning and end of every iteration is equal to the
     *    total violation of the sum condition.
     */
    int count = 0;
    while (true) {
        /*
         * run DFS to compute a pre-order traversal of the tree
         * to enable constant time ancestor queries
         */
        std::stack<int> call_stack;
        call_stack.push(root);
        std::vector<int> order(clone_tree.nodes().size());
        int order_index = 0;
        while(!call_stack.empty()) {
            int i = call_stack.top();
            call_stack.pop();

            order[i] = order_index;
            order_index++;

            for (auto k : clone_tree.successors(i)) {
                call_stack.push(k);
            }
        }

        /* 
         * Construct a priority queue of all possible moves, ranked by 
         * the amount they reduce the total violation score.  Since the 
         * total violation is a lower bound on the reconstruction error,
         * we can stop when the total violation is greater than the 
         * best score.
         */
        std::set<std::tuple<float, int, int>> move_priority_queue;
        for (auto u : clone_tree.nodes()) { // quadratic iteration complexity: okay?
            for (auto v : clone_tree.nodes()) {
                if (u == v || u == root || order[u] <= order[v]) {
                    continue; 
                }

                int p = *clone_tree.predecessors(u).begin();
                int p_column = clone_tree[p].data.id;
                int u_column = clone_tree[u].data.id;
                int v_column = clone_tree[v].data.id;

                float diff = 0;
                for (size_t j = 0; j < nrows; ++j) {
                    float A_jp = A[p_column][j] - F_transpose[u_column][j];
                    diff += std::max(A_jp - F_transpose[p_column][j], 0.0f) - std::max(A[p_column][j] - F_transpose[p_column][j], 0.0f);
                }

                for (size_t j = 0; j < nrows; ++j) {
                    float A_jv = A[v_column][j] + F_transpose[u_column][j];
                    diff += std::max(A_jv - F_transpose[v_column][j], 0.0f) - std::max(A[v_column][j] - F_transpose[v_column][j], 0.0f);
                }

                move_priority_queue.insert({total_violation + diff, u, v});
           }
        }

        float best_score = one_fastbe(clone_tree, vertex_map, F, root);

        if (count % progress_interval == 0) {
            spdlog::info("Sample ID {}: iteration {}, objective value: {}, total violation: {}", sample_id, count, best_score, total_violation);
        }

        std::pair<int, int> best_move = {-1, -1};
        int iterations = 0;
        for (auto [tv, u, v] : move_priority_queue) {
            if (tv > best_score || iterations >= max_iterations) {
                break;
            }

            iterations++;
        
            int parent = *clone_tree.predecessors(u).begin();
            subtree_prune_and_regraft(clone_tree, u, v, root);

            float score = one_fastbe(clone_tree, vertex_map, F, root);
            subtree_prune_and_regraft(clone_tree, u, parent, root);

            if (score < best_score - 1e-4) {
                best_score = score;
                best_move = {u, v};
            } 
        }

        if (best_move.first == -1) {
            spdlog::info("Sample ID {} found local minimum -- no SPR moves can improve the score.", sample_id);
            break;
        }

        auto [u, v] = best_move;
        subtree_prune_and_regraft(clone_tree, u, v, root);

        A = compute_sum_violation_matrix(clone_tree, vertex_map, F_transpose);
        total_violation = 0;
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < A[0].size(); ++j) {
                total_violation += std::max(A[i][j] - F_transpose[i][j], 0.0f);
            }
        }

        count++;
    }

    spdlog::info("Sample ID {} hill climbing completed @ {} iterations", sample_id, count);
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

void perform_regression(argparse::ArgumentParser regress) {
    auto [clone_tree, vertex_map] = parse_adjacency_list(regress.get<std::string>("clone_tree"));
    std::vector<std::vector<float>> frequency_matrix = parse_frequency_matrix(regress.get<std::string>("frequency_matrix"));

    std::vector<int> columns(vertex_map.size());
    for (auto [v, i] : vertex_map) {
        columns[i] = v;
    }

    auto start = std::chrono::high_resolution_clock::now();
    float obj = one_fastbe(clone_tree, vertex_map, frequency_matrix, columns, 0);
    for (size_t i = 0; i < regress.get<size_t>("num_reps") - 1; ++i) {
        one_fastbe(clone_tree, vertex_map, frequency_matrix, columns, 0);
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

void perform_search(argparse::ArgumentParser search) {
    std::vector<std::vector<float>> frequency_matrix = parse_frequency_matrix(search.get<std::string>("frequency_matrix"));

    size_t nrows = frequency_matrix.size();
    size_t ncols = frequency_matrix[0].size();

    float ancestry_threshold = search.get<float>("ancestry_threshold");
    digraph<int> ancestry_graph;
    for (size_t i = 0; i < ncols; ++i) {
        ancestry_graph.add_vertex(i);
    }

    int assigned_root = search.get<int>("assigned_root");

    size_t num_edges = 0;
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            if (i == j) continue;

            bool to_add = true;
            for (size_t k = 0; k < nrows; ++k) {
                // violates ancestry condition if the frequency of i is less than the frequency of j
                if (frequency_matrix[k][i] <= frequency_matrix[k][j] - ancestry_threshold) {
                    to_add = false;
                }
            }

            if (to_add) {
                ancestry_graph.add_edge(i, j);
                ++num_edges;
            }
        }
    }

    std::map<std::pair<int,int>, double> G_weights;
    for (auto [i, j] : ancestry_graph.edges()) {
        G_weights[std::make_pair(i, j)] = 1.0;
    }

    spdlog::info("Number of edges in ancestry graph: {}", num_edges);

    std::unordered_map<int, int> identity_map;
    for (size_t i = 0; i < ncols; ++i) {
        identity_map[i] = i;
    }

    if (search.get<bool>("benchmark")) {
        spdlog::info("Benchmarking computing the l1 minimizer 100 times...");
        spdlog::info("Number of samples (m={}) and  number of clones (n={})", nrows, ncols);


        // use fixed seed for benchmarking
        // std::random_device rd;
        // std::ranlux48_base gen(rd());
        std::ranlux48_base gen(0);

        auto [clone_tree_int, root] = sample_random_spanning_tree(ancestry_graph, G_weights, gen, assigned_root);
        digraph<clone_tree_vertex> clone_tree = convert_clone_tree(clone_tree_int, nrows);
        auto start = std::chrono::high_resolution_clock::now();
        float obj1 = one_fastbe(clone_tree, identity_map, frequency_matrix, root);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        std::cout << "Time to compute l1 minimizer: " << duration.count() << " ns" << std::endl;

        spdlog::info("Objective value (one_fastbe): {}", obj1);

        return;
    }

    /* Draw N random spanning trees from the complete graph and then
     * perform hill climbing on each of them. */
    std::mutex mtx;
    std::atomic<size_t> counter(0);

    size_t num_samples = search.get<size_t>("samples");
    unsigned int num_threads = search.get<unsigned int>("threads");
    std::vector<std::thread> threads(num_threads);

    std::pair<digraph<clone_tree_vertex>, float> best_tree = std::make_pair(digraph<clone_tree_vertex>(), std::numeric_limits<float>::max());
    for (auto& th : threads) {
        th = std::thread([&]() {
            size_t thread_id = std::hash<std::thread::id>{}(std::this_thread::get_id());
            (void) thread_id; // suppress unused variable warning
            while (true) {
                std::random_device rd;
                std::ranlux48_base gen(rd());
                size_t i = counter++;
                if (i >= num_samples) return;  // no more work to do

                /* Create a random tree by picking a random spanning tree of G. */
                auto [clone_tree_int, root] = sample_random_spanning_tree(ancestry_graph, G_weights, gen, assigned_root);
                digraph<clone_tree_vertex> clone_tree = convert_clone_tree(clone_tree_int, nrows);

                if (search.get<bool>("two-phase")) {
                    spdlog::info("Sample ID {}: performing total violation hill climb (phase I)...", i);
                    deterministic_total_violation_hill_climb(clone_tree, identity_map, frequency_matrix, root, ncols, i, search.get<int>("progress_interval"));

                    spdlog::info("Sample ID {}: performing hill climb (phase II)...", i);
                    deterministic_hill_climb(clone_tree, identity_map, frequency_matrix, root, ncols, i, search.get<int>("progress_interval"));
                } else {
                    spdlog::info("Sample ID {}: performing hill climb...", i);
                    deterministic_hill_climb(clone_tree, identity_map, frequency_matrix, root, ncols, i, search.get<int>("progress_interval"));
                }

                float obj = one_fastbe(clone_tree, identity_map, frequency_matrix, root);

                spdlog::info("Random hill climbed tree objective: {}", obj);

                {
                    std::lock_guard<std::mutex> lock(mtx);  // lock the mutex to protect shared resources
                    spdlog::info("Random sample iteration {}, Objective value: {}", i, obj);
                    if (obj < best_tree.second) {
                        best_tree = std::make_pair(clone_tree, obj);
                    }
                }
            }
        });
    }

    for (auto& th : threads) {
        th.join();  // wait for all threads to finish
    }

    digraph<clone_tree_vertex> min_tree = best_tree.first;
    float min_obj = best_tree.second;

    std::string adjacency_list = to_adjacency_list(min_tree, identity_map);

    std::ofstream adj_output(search.get<std::string>("output") + "_tree.txt");
    adj_output << adjacency_list;

    json output;
    output["objective_value"] = min_obj;
    std::ofstream json_output(search.get<std::string>("output") + "_results.json");
    json_output << output.dump(4) << std::endl;

    return;
}

int main(int argc, char *argv[])
{
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

    regress.add_description("regresses a clone tree onto a frequency matrix.");
    search.add_description("searches for a clone tree that best fits a frequency matrix.");

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

    search.add_argument("frequency_matrix")
          .help("TXT file containing the frequency matrix");

    search.add_argument("-o", "--output")
          .help("prefix of the output files")
          .required();

    search.add_argument("-s", "--samples")
          .help("number of random samples to use for hill climbing")
          .default_value((size_t) 100)
          .scan<'u', size_t>();

    search.add_argument("-r", "--random-seed") // TODO: get working
          .help("random seed")
          .default_value(0)
          .scan<'d', int>();

    search.add_argument("-t", "--threads")
          .help("number of threads to use")
          .default_value(std::thread::hardware_concurrency())
          .scan<'u', unsigned int>();

    search.add_argument("-a", "--ancestry_threshold")
          .help("ancestry threshold")
          .default_value(1.0f)
          .scan<'g', float>();

    search.add_argument("-f", "--assigned_root")
          .help("assigned root")
          .default_value(0)
          .scan<'d', int>();

    search.add_argument("-p", "--progress_interval")
          .help("progress interval")
          .default_value(5)
          .scan<'d', int>();

    search.add_argument("--benchmark")
          .help("benchmark mode")
          .default_value(false)
          .implicit_value(true);

    search.add_argument("--two-phase")
          .help("two phase mode")
          .default_value(false)
          .implicit_value(true);

    program.add_subparser(search);
    program.add_subparser(regress);
    
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;

        if (program.is_subcommand_used(regress)) {
            std::cerr << regress;
        } else if (program.is_subcommand_used(search)) {
            std::cerr << search;
        } else {
            std::cerr << program;
        }

        std::exit(1);
    }

    if (program.is_subcommand_used(regress)) {
        perform_regression(regress);
    } else if (program.is_subcommand_used(search)) {
        perform_search(search);
    } else {
        std::cerr << program;
    }

    return 0;
}
