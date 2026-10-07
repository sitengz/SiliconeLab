#ifndef SILICONE_CONFIG_INPUT_HPP
#define SILICONE_CONFIG_INPUT_HPP

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace silicone_config {

inline std::string trim(const std::string& value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

inline std::string option_name(std::string key) {
    key = trim(key);
    if (key.rfind("--", 0) == 0) key.erase(0, 2);
    std::replace(key.begin(), key.end(), '_', '-');
    return "--" + key;
}

template <typename Apply>
void read_config_file(const std::string& path, Apply apply) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open config file: " + path);

    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        line = trim(line);
        if (line.empty()) continue;

        const std::size_t separator = line.find('=');
        if (separator == std::string::npos)
            throw std::runtime_error(
                "Config line " + std::to_string(line_number) +
                " must use key = value syntax");
        const std::string key = option_name(line.substr(0, separator));
        std::string value = trim(line.substr(separator + 1));
        if (value.size() >= 2 &&
            ((value.front() == '"' && value.back() == '"') ||
             (value.front() == '\'' && value.back() == '\'')))
            value = value.substr(1, value.size() - 2);
        if (value.empty())
            throw std::runtime_error(
                "Config line " + std::to_string(line_number) +
                " has an empty value");
        try {
            apply(key, value);
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Config line " + std::to_string(line_number) + ": " +
                error.what());
        }
    }
    if (!input.good() && !input.eof())
        throw std::runtime_error("Failed while reading config file: " + path);
}

template <typename Apply, typename Help>
std::string parse_arguments(int argc, char** argv, Apply apply, Help help) {
    std::vector<std::pair<std::string, std::string>> options;
    std::string config_path;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help") {
            help();
            std::exit(0);
        }
        if (i + 1 >= argc)
            throw std::runtime_error("Missing value after " + option);
        const std::string value = argv[++i];
        if (option == "--config") {
            if (!config_path.empty())
                throw std::runtime_error("--config may be supplied only once");
            config_path = value;
        } else {
            options.emplace_back(option, value);
        }
    }

    if (!config_path.empty()) read_config_file(config_path, apply);
    for (const auto& option : options) apply(option.first, option.second);
    return config_path;
}

} // namespace silicone_config

#endif
