#include <pprint.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>

#include <vafpp.hpp>
#include <digraph.hpp>

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
void perform_regression(argparse::ArgumentParser regress) {
    digraph<int> clone_tree = parse_adjacency_list(regress.get<std::string>("clone_tree"));
    for (auto v : clone_tree.nodes()) {
        spdlog::info("Vertex {} has {} children", clone_tree[v].data, clone_tree.out_degree(v));
        for (auto child : clone_tree.successors(v)) {
            spdlog::info("Child {}", clone_tree[child].data);
        }
    }

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
