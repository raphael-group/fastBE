#include <pprint.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <argparse/argparse.hpp>

#include <vafpp.hpp>
#include <digraph.hpp>
#include <piecewiselinearf.hpp>

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

using namespace std;
using json = nlohmann::json;

/*
 * Given a frequency matrix $F$ and a clone tree $T$, this function
 * finds the minimizing value of $$\sum_{i=1}^m\lVert F_i - (UB)_i \rVert_1$$ 
 * over all usage matrices in $\mathcal{O}(mn^2)$ using dynamic programming 
 * and an analysis of convex and continuous,  piecewise linear functions.
*/
double one_vafpp(const digraph<int> clone_tree, std::vector<std::vector<double>>& F) {
    size_t nrows = F.size();
    size_t ncols = F[0].size();

    std::vector<std::vector<double>> W(nrows, std::vector<double>(ncols, 0.0));

    for (size_t j = 0; j < nrows; ++j) {
        for (size_t i = 0; i < ncols; ++i) {
            W[j][i] = F[j][i];
            for (auto k : clone_tree.successors(i)) {
                W[j][i] -= F[j][clone_tree[k].data];
            }
        }
    }

    std::function<PiecewiseLinearF(size_t, size_t)> one_vafpp_recursive = [&](size_t j, size_t i) {
        if (clone_tree.out_degree(i) == 0) {
            return PiecewiseLinearF({W[j][i]}, 0); 
        }

        PiecewiseLinearF g_out({W[j][i]}, 0);
        for (auto k : clone_tree.successors(i)) {
            auto f = one_vafpp_recursive(j, clone_tree[k].data);
            auto g = compute_minimizer(f);
            g_out = g_out + g;
        }

        return g_out;
    };

    double obj = 0;
    for (size_t j = 0; j < nrows; ++j) {
        auto f = one_vafpp_recursive(j, 0);
        f = compute_minimizer(f) + PiecewiseLinearF({1 - F[j][0]}, 0);
        obj += *min_element(f.intercepts.begin(), f.intercepts.end());
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
    digraph<int> clone_tree = parse_adjacency_list(regress.get<std::string>("clone_tree"));
    std::vector<std::vector<double>> frequency_matrix = parse_frequency_matrix(regress.get<std::string>("frequency_matrix"));

    double obj = one_vafpp(clone_tree, frequency_matrix);
    std::cout << obj << std::endl;
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

    regress.add_description("Regresses a clone tree onto a frequency matrix.");

    regress.add_argument("clone_tree")
           .help("adjacency list of the clone tree");

    regress.add_argument("frequency_matrix")
           .help("TXT file containing the frequency matrix");

    regress.add_argument("-o", "--output")
           .help("prefix of the output files")
           .required();

    regress.add_argument("-a", "--aggression")
           .help("aggression of stochastic perturbation in (0, infinity)")
           .default_value(1.0)
           .scan<'g', double>();

    regress.add_argument("-i", "--iterations")
           .help("number of iterations to perform without improvement before stopping")
           .default_value(100)
           .scan<'d', int>();

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
    } else {
        std::cerr << program;
    }

    return 0;
}
