#include <pprint.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>
#include <vafpp.hpp>

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

void do_search(argparse::ArgumentParser search) {
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

    argparse::ArgumentParser search(
        "search"
    );

    search.add_description("Infers a clone tree using local tree rearrangement operations");

    search.add_argument("clone_tree")
          .help("adjacency list of the clone tree");

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
          .scan<'d', int>();

    program.add_subparser(search);
    
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;

        if (program.is_subcommand_used(search)) {
            std::cerr << search;
        } else {
            std::cerr << program;
        }

        std::exit(1);
    }

    if (program.is_subcommand_used(search)) {
        do_search(search);
    } else {
        std::cerr << program;
    }

    return 0;
}
