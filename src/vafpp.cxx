#include <pprint.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <argparse/argparse.hpp>

#include <vafpp.hpp>
#include <digraph.hpp>
#include <piecewiselinearf.hpp>

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
 * Given a frequency matrix $F$ and a clone tree $T$, this function
 * finds the minimizing value of $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$ 
 * over all usage matrices in $\mathcal{O}(mn^2)$ using dynamic programming 
 * and an analysis of convex and continuous, piecewise linear functions.
 *
 * Input: 
 *  - clone_tree: A clone tree represented as a digraph.
 *  - vertex_map: A map from the rows of the frequency matrix to the 
 *    vertices of the clone tree.
 *  - F: A frequency matrix represented as a 2D vector.
*/
double one_vafpp(const digraph<int>& clone_tree, const std::map<int, int>& vertex_map, std::vector<std::vector<double>>& F, int root) {
    size_t nrows = F.size();
    size_t ncols = F[0].size();

    std::vector<std::vector<double>> W(nrows, std::vector<double>(ncols, 0.0));

    for (size_t j = 0; j < nrows; ++j) {
        for (size_t i = 0; i < ncols; ++i) {
            W[j][i] = F[j][i];
            for (auto k : clone_tree.successors(vertex_map.at(i))) {
                W[j][i] -= F[j][clone_tree[k].data];
            }
        }
    }

    std::function<PiecewiseLinearF(size_t, size_t)> one_vafpp_recursive = [&](size_t j, size_t i) {
        if (clone_tree.out_degree(vertex_map.at(i)) == 0) {
            return PiecewiseLinearF({W[j][i]}, 0); 
        }

        PiecewiseLinearF g_out({W[j][i]}, 0);
        for (auto k : clone_tree.successors(vertex_map.at(i))) {
            auto f = one_vafpp_recursive(j, clone_tree[k].data);
            f.compute_minimizer();
            g_out.addInPlace(std::move(f));
        }

        return g_out;
    };

    double obj = 0;
    for (size_t j = 0; j < nrows; ++j) {
        auto f = one_vafpp_recursive(j, root);
        f.compute_minimizer();
        f.addInPlace(PiecewiseLinearF({1 - F[j][0]}, 0));
        obj += f.minimizer();
    }

    return -1 * obj;
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
std::vector<std::vector<double>> parse_frequency_matrix(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open the file.");
    }

    std::vector<std::vector<double>> matrix;
    std::string line;
    size_t num_cols = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        double value;

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
    std::vector<std::vector<double>> frequency_matrix = parse_frequency_matrix(regress.get<std::string>("frequency_matrix"));

    auto start = std::chrono::high_resolution_clock::now();
    double obj = one_vafpp(clone_tree, vertex_map, frequency_matrix, 0);

    for (size_t i = 0; i < regress.get<size_t>("num_reps") - 1; ++i) {
        one_vafpp(clone_tree, vertex_map, frequency_matrix, 0);
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
    std::vector<std::vector<double>> frequency_matrix = parse_frequency_matrix(search.get<std::string>("frequency_matrix"));

    size_t nrows = frequency_matrix.size();
    size_t ncols = frequency_matrix[0].size();

    digraph<int> complete_graph;
    for (size_t i = 0; i < ncols; ++i) {
        complete_graph.add_vertex(i);
    }

    for (size_t i = 0; i < ncols; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            if (i == j) continue;
            complete_graph.add_edge(i, j);
        }
    }


    std::random_device rd;
    std::ranlux48_base gen(rd());

    std::map<int, int> identity_map;
    for (size_t i = 0; i < ncols; ++i) {
        identity_map[i] = i;
    }

    float min_obj = std::numeric_limits<float>::infinity();
    std::optional<digraph<int>> min_tree;
    for (size_t i = 0; i < search.get<size_t>("iterations"); i++) {
        auto [clone_tree, root] = sample_random_spanning_tree(complete_graph, gen);
        float obj = one_vafpp(clone_tree, identity_map, frequency_matrix, root);
        if (obj < min_obj) {
            min_obj = obj;
            min_tree = clone_tree;
        }
    }

    std::string adjacency_list = to_adjacency_list(*min_tree, identity_map);

    std::ofstream adj_output(search.get<std::string>("output") + "_tree.txt");
    adj_output << adjacency_list << std::endl;

    json output;
    output["objective_value"] = min_obj;
    std::ofstream json_output(search.get<std::string>("output") + "_results.json");
    json_output << output.dump(4) << std::endl;

    return;
}

int main(int argc, char *argv[])
{
    auto console_logger = spdlog::stdout_color_mt("vafpp");
    spdlog::set_default_logger(console_logger);

    auto error_logger = spdlog::stderr_color_mt("error");

    argparse::ArgumentParser program(
        "vafpp",
        std::to_string(VAFPP_VERSION_MAJOR) + "." + std::to_string(VAFPP_VERSION_MINOR)
    );

    argparse::ArgumentParser regress(
        "regress"
    );

    argparse::ArgumentParser search(
        "search"
    );

    regress.add_description("Regresses a clone tree onto a frequency matrix.");

    regress.add_argument("clone_tree")
           .help("adjacency list of the clone tree");

    regress.add_argument("frequency_matrix")
           .help("TXT file containing the frequency matrix");

    regress.add_argument("-o", "--output")
           .help("prefix of the output files")
           .required();

    regress.add_argument("-n", "--num_reps")
           .help("number of times to repeat the regression for benchmarking")
           .default_value(1)
           .scan<'d', size_t>();

    search.add_argument("frequency_matrix")
          .help("TXT file containing the frequency matrix");

    search.add_argument("-o", "--output")
          .help("prefix of the output files")
          .required();

    search.add_argument("-a", "--aggression")
          .help("aggression of stochastic perturbation in (0, infinity)")
          .default_value(1.0)
          .scan<'g', double>();

    search.add_argument("-i", "--iterations")
          .help("number of iterations to perform without improvement before stopping")
          .default_value(100)
          .scan<'d', size_t>();

    program.add_subparser(search);
    program.add_subparser(regress);
    
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;

        if (program.is_subcommand_used(regress)) {
            std::cerr << regress;
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
