#include "utils.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <filesystem>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <path>\n";
        return 1;
    }

    const std::string path(argv[1]);

    // Extract m and n from path component like "m32_n32_k32"
    std::regex dim_re(R"(m(\d+)_n(\d+))");
    std::smatch match;
    if (!std::regex_search(path, match, dim_re)) {
        std::cerr << "Error: could not extract m and n from path: " << path << "\n";
        return 1;
    }

    const int m = std::stoi(match[1]);
    const int n = std::stoi(match[2]);

    std::cout << "m=" << m << " n=" << n << "\n";

    std::vector<double> data(m * n);
    if (!load_matrix(path, m, n, &data[0], m)) {
        std::cerr << "Error: load_matrix failed for: " << path << "\n";
        return 1;
    }

    display_distribution_stats(m, n, &data[0], m, path);

    return 0;
}