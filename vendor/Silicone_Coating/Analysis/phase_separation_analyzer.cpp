#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

enum ComponentIndex {
    kNetwork = 0,
    kCrosslinker = 1,
    kFiller = 2,
    kModerator = 3,
    kComponentCount = 4
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Box {
    std::array<double, 3> lo{{0.0, 0.0, 0.0}};
    std::array<double, 3> hi{{0.0, 0.0, 0.0}};
    std::array<bool, 3> periodic{{false, false, false}};

    double length(std::size_t dimension) const {
        return hi[dimension] - lo[dimension];
    }
};

struct ComponentInfo {
    std::string key;
    long long molecules = 0;
    long long beads = 0;
};

struct ModelInfo {
    std::string format;
    long long format_version = 0;
    std::string formulation;
    std::string case_name;
    std::string geometry;
    std::array<ComponentInfo, kComponentCount> components;
    long long total_beads = 0;
    long long total_molecules = 0;
    std::string oil_model;
    long long dms_repeats_per_chain = 0;
    long long mps_repeats_per_chain = 0;
    long long total_mps_repeats = 0;
    long long chains_with_extra_mps = 0;
    long long oil_chain_count = 0;
    double timestep_fs = 0.0;
    long long production_steps = 0;
    long long dump_every_steps = 0;
    long long declared_expected_frames = 0;
    long long expected_frames = 0;
};

struct Options {
    std::string input_file;
    std::string info_file;
    long long frame_stride = 10;
    double target_grid_spacing = 8.0;
    double contact_cutoff = 8.0;
    double q_max = std::numeric_limits<double>::quiet_NaN();
    std::string metrics_file;
    std::string structure_file;
    std::string field_file;
    std::string report_file;
};

struct Marker {
    long long atom = 0;
    long long molecule = 0;
    Vec3 position;
};

struct Frame {
    long long index = 0;
    long long timestep = 0;
    Box box;
    std::vector<Marker> dms;
    std::vector<Marker> mps;
    long long pendant_count = 0;
};

struct CompositionGrid {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    Box box;
    std::vector<long long> mps;
    std::vector<long long> dms;

    std::size_t index(int ix, int iy, int iz) const {
        return (static_cast<std::size_t>(iz) * static_cast<std::size_t>(ny) +
                static_cast<std::size_t>(iy)) *
                   static_cast<std::size_t>(nx) +
               static_cast<std::size_t>(ix);
    }
};

struct SegregationResult {
    double observed_variance = 0.0;
    double random_variance = 0.0;
    double index = 0.0;
    long long occupied_cells = 0;
};

struct SpectrumBin {
    double q = 0.0;
    long long modes = 0;
    double scc = 0.0;
};

struct SpectrumResult {
    std::vector<SpectrumBin> bins;
    double q_max_used = 0.0;
    double dq = 0.0;
    double low_q_mean = 0.0;
    double q_peak = 0.0;
    double peak_scc = 0.0;
    double domain_size = 0.0;
    bool peak_at_lowest_q = false;
};

struct ClusterResult {
    long long clusters = 0;
    long long largest_cluster_chains = 0;
    double largest_cluster_fraction = 0.0;
    long long unique_chain_contacts = 0;
};

struct PhaseMetrics {
    long long frame_index = 0;
    long long timestep = 0;
    double time_ns = 0.0;
    long long mps_markers = 0;
    long long dms_markers = 0;
    double mps_fraction = 0.0;
    SegregationResult segregation_3d;
    SegregationResult segregation_xy;
    SpectrumResult spectrum;
    ClusterResult clusters;
};

std::string trim(const std::string &input) {
    const std::string whitespace = " \t\r\n";
    const std::size_t first = input.find_first_not_of(whitespace);
    if (first == std::string::npos) return "";
    const std::size_t last = input.find_last_not_of(whitespace);
    return input.substr(first, last - first + 1);
}

std::vector<std::string> split_words(const std::string &line) {
    std::istringstream input(line);
    std::vector<std::string> words;
    std::string word;
    while (input >> word) words.push_back(word);
    return words;
}

std::string read_text_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open info file: " + path);
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading info file: " + path);
    }
    return contents.str();
}

std::size_t matching_brace(const std::string &text, std::size_t open) {
    if (open >= text.size() || text[open] != '{') {
        throw std::runtime_error("internal JSON parsing error");
    }
    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t index = open; index < text.size(); ++index) {
        const char value = text[index];
        if (in_string) {
            if (escaped) escaped = false;
            else if (value == '\\') escaped = true;
            else if (value == '"') in_string = false;
            continue;
        }
        if (value == '"') in_string = true;
        else if (value == '{') ++depth;
        else if (value == '}' && --depth == 0) return index;
    }
    throw std::runtime_error("unterminated JSON object");
}

std::string json_object_for_key(const std::string &text, const std::string &key) {
    const std::string quoted = "\"" + key + "\"";
    const std::size_t key_position = text.find(quoted);
    if (key_position == std::string::npos) {
        throw std::runtime_error("missing JSON object: " + key);
    }
    const std::size_t colon = text.find(':', key_position + quoted.size());
    const std::size_t open = colon == std::string::npos
        ? std::string::npos : text.find('{', colon + 1);
    if (open == std::string::npos) {
        throw std::runtime_error("JSON field is not an object: " + key);
    }
    const std::size_t close = matching_brace(text, open);
    return text.substr(open, close - open + 1);
}

std::vector<std::pair<std::string, std::string>>
json_child_objects(const std::string &object_text) {
    std::vector<std::pair<std::string, std::string>> children;
    std::size_t position = 1;
    while (position + 1 < object_text.size()) {
        while (position < object_text.size() &&
               (std::isspace(static_cast<unsigned char>(object_text[position])) ||
                object_text[position] == ',')) {
            ++position;
        }
        if (position >= object_text.size() || object_text[position] == '}') break;
        if (object_text[position] != '"') {
            throw std::runtime_error("invalid key in JSON object");
        }
        const std::size_t key_begin = ++position;
        while (position < object_text.size() && object_text[position] != '"') {
            if (object_text[position] == '\\') ++position;
            ++position;
        }
        if (position >= object_text.size()) {
            throw std::runtime_error("unterminated JSON key");
        }
        const std::string key = object_text.substr(key_begin, position - key_begin);
        const std::size_t colon = object_text.find(':', ++position);
        const std::size_t open = colon == std::string::npos
            ? std::string::npos
            : object_text.find_first_not_of(" \t\r\n", colon + 1);
        if (open == std::string::npos || object_text[open] != '{') {
            throw std::runtime_error("JSON child is not an object: " + key);
        }
        const std::size_t close = matching_brace(object_text, open);
        children.emplace_back(
            key, object_text.substr(open, close - open + 1));
        position = close + 1;
    }
    return children;
}

double json_number(const std::string &text, const std::string &key) {
    const std::regex pattern(
        "\"" + key +
        "\"\\s*:\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        throw std::runtime_error("missing numeric info field: " + key);
    }
    return std::stod(match[1].str());
}

long long json_integer(const std::string &text, const std::string &key) {
    const double value = json_number(text, key);
    if (!std::isfinite(value) || value < 0.0 ||
        std::fabs(value - std::round(value)) > 1.0e-8) {
        throw std::runtime_error(
            "info field is not a nonnegative integer: " + key);
    }
    return static_cast<long long>(std::llround(value));
}

std::string json_string(const std::string &text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        throw std::runtime_error("missing string info field: " + key);
    }
    return match[1].str();
}

ComponentInfo component_from_entry(
    const std::pair<std::string, std::string> &entry) {
    ComponentInfo component;
    component.key = entry.first;
    component.molecules = json_integer(entry.second, "M");
    component.beads = json_integer(entry.second, "beads");
    return component;
}

ModelInfo parse_model_info(const std::string &path) {
    const std::string text = read_text_file(path);
    ModelInfo info;
    info.format = json_string(text, "format");
    info.format_version = json_integer(text, "format_version");
    info.formulation = json_string(text, "formulation");
    info.case_name = json_string(text, "case_name");
    info.geometry = json_string(text, "geometry");
    if (info.format_version != 2 ||
        (info.format != "V22-model-info" && info.format != "V35-model-info")) {
        throw std::runtime_error(
            "phase analysis requires a current V22/V35 model-info version-2 file");
    }
    if (info.geometry != "bulk" && info.geometry != "film") {
        throw std::runtime_error("unsupported geometry: " + info.geometry);
    }

    const auto component_entries =
        json_child_objects(json_object_for_key(text, "components"));
    bool have_network = false;
    bool have_crosslinker = false;
    bool have_filler = false;
    bool have_moderator = false;
    for (const auto &entry : component_entries) {
        if (entry.first == "network_strands") {
            info.components[kNetwork] = component_from_entry(entry);
            have_network = true;
        } else if (entry.first == "crosslinkers") {
            info.components[kCrosslinker] = component_from_entry(entry);
            have_crosslinker = true;
        } else if (entry.first == "star_moderators") {
            info.components[kModerator] = component_from_entry(entry);
            have_moderator = true;
        } else {
            if (have_filler) {
                throw std::runtime_error(
                    "more than one filler component in info file");
            }
            info.components[kFiller] = component_from_entry(entry);
            have_filler = true;
        }
    }
    if (!have_network || !have_crosslinker || !have_filler || !have_moderator) {
        throw std::runtime_error(
            "could not identify all four formulation components");
    }

    const std::string composition = json_object_for_key(text, "composition");
    info.total_beads = json_integer(composition, "total_beads");
    info.total_molecules = json_integer(composition, "total_molecules");
    long long bead_sum = 0;
    long long molecule_sum = 0;
    for (const ComponentInfo &component : info.components) {
        bead_sum += component.beads;
        molecule_sum += component.molecules;
    }
    if (bead_sum != info.total_beads ||
        molecule_sum != info.total_molecules) {
        throw std::runtime_error(
            "component totals do not match composition totals");
    }

    const std::string oil = json_object_for_key(text, "silicone_oil");
    info.oil_model = json_string(oil, "model");
    info.dms_repeats_per_chain =
        json_integer(oil, "dms_repeats_per_chain");
    info.mps_repeats_per_chain =
        json_integer(oil, "mps_repeats_per_chain");
    info.oil_chain_count = json_integer(oil, "chain_count");
    info.total_mps_repeats = oil.find("\"total_mps_repeats\"") != std::string::npos
        ? json_integer(oil, "total_mps_repeats")
        : info.oil_chain_count * info.mps_repeats_per_chain;
    if (oil.find("\"chains_with_extra_mps\"") != std::string::npos)
        info.chains_with_extra_mps = json_integer(oil, "chains_with_extra_mps");
    if (info.total_mps_repeats !=
        info.oil_chain_count * info.mps_repeats_per_chain +
        info.chains_with_extra_mps)
        throw std::runtime_error("silicone-oil MPS totals are inconsistent");
    if (info.oil_chain_count != info.components[kFiller].molecules) {
        throw std::runtime_error(
            "silicone-oil chain count does not match filler component");
    }

    const std::string simulation =
        json_object_for_key(text, "simulation_template");
    info.timestep_fs = json_number(simulation, "timestep_fs");
    const std::string production =
        json_object_for_key(simulation, "msd_production");
    info.production_steps = json_integer(production, "steps");
    info.dump_every_steps = json_integer(production, "dump_every_steps");
    info.declared_expected_frames =
        json_integer(production, "expected_frames");
    if (info.timestep_fs <= 0.0 || info.dump_every_steps <= 0) {
        throw std::runtime_error("invalid trajectory timing fields in info file");
    }
    info.expected_frames =
        info.production_steps / info.dump_every_steps + 1;
    return info;
}

std::array<long long, kComponentCount> component_ends(const ModelInfo &info) {
    std::array<long long, kComponentCount> ends{{0, 0, 0, 0}};
    long long total = 0;
    for (std::size_t component = 0; component < kComponentCount; ++component) {
        total += info.components[component].molecules;
        ends[component] = total;
    }
    return ends;
}

int component_for_molecule(
    long long molecule,
    const std::array<long long, kComponentCount> &ends) {
    if (molecule < 1 || molecule > ends.back()) {
        throw std::runtime_error(
            "trajectory molecule ID outside info ranges: " +
            std::to_string(molecule));
    }
    for (int component = 0; component < kComponentCount; ++component) {
        if (molecule <= ends[static_cast<std::size_t>(component)]) {
            return component;
        }
    }
    throw std::runtime_error("internal molecule-range error");
}

std::string sanitize_name(std::string name) {
    for (char &value : name) {
        const unsigned char byte = static_cast<unsigned char>(value);
        if (!std::isalnum(byte) && value != '-' && value != '_' && value != '.') {
            value = '_';
        }
    }
    if (name.empty()) name = "sample";
    return name;
}

void resolve_output_paths(Options &options, const ModelInfo &info) {
    namespace fs = std::filesystem;
    const std::string sample = sanitize_name(info.case_name);
    const fs::path folder(sample);
    if (options.metrics_file.empty()) {
        options.metrics_file =
            (folder / ("phase_metrics." + sample + ".dat")).string();
    }
    if (options.structure_file.empty()) {
        options.structure_file =
            (folder / ("phase_structure_factor." + sample + ".dat")).string();
    }
    if (options.field_file.empty()) {
        options.field_file =
            (folder / ("phase_field_final." + sample + ".dat")).string();
    }
    if (options.report_file.empty()) {
        options.report_file =
            (folder / ("phase_report." + sample + ".txt")).string();
    }
    const std::array<std::string, 4> paths{{
        options.metrics_file,
        options.structure_file,
        options.field_file,
        options.report_file
    }};
    for (const std::string &path : paths) {
        const fs::path parent = fs::path(path).parent_path();
        if (!parent.empty()) fs::create_directories(parent);
    }
}

long long parse_positive_integer(
    const std::string &value,
    const std::string &option) {
    std::size_t consumed = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(value, &consumed);
    } catch (...) {
        throw std::runtime_error(option + " requires a positive integer");
    }
    if (consumed != value.size() || parsed <= 0) {
        throw std::runtime_error(option + " requires a positive integer");
    }
    return parsed;
}

double parse_positive_double(
    const std::string &value,
    const std::string &option) {
    std::size_t consumed = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (...) {
        throw std::runtime_error(option + " requires a positive number");
    }
    if (consumed != value.size() || !std::isfinite(parsed) || parsed <= 0.0) {
        throw std::runtime_error(option + " requires a positive number");
    }
    return parsed;
}

void print_help(const char *program) {
    std::cout
        << "Usage:\n  " << program << " DATA_OR_TRAJECTORY INFO_FILE [options]\n\n"
        << "Measure DMS/MPS phase separation in a final V22/V35 data snapshot\n"
        << "or the generated MSD trajectory.\n\n"
        << "Options:\n"
        << "  --frame-stride N       analyze every Nth frame; default: 10\n"
        << "  --grid-spacing X       target composition-grid spacing in A; default: 8\n"
        << "  --contact-cutoff X     intermolecular MPS contact cutoff in A; default: 8\n"
        << "  --q-max X              largest structure-factor q in 1/A; default: auto\n"
        << "  --output FILE          phase-metrics table path\n"
        << "  --structure-output FILE  concentration structure-factor table\n"
        << "  --field-output FILE    final 3D composition-field table\n"
        << "  --report-output FILE   text report path\n"
        << "  --help                 show this reference\n\n"
        << "Input type is detected automatically. For a trajectory, the final\n"
        << "frame is always analyzed even when it is not on the requested stride.\n"
        << "Bulk systems use 3D S_cc(q); films use in-plane S_cc(q_xy).\n";
}

Options parse_options(int argc, char **argv) {
    Options options;
    std::vector<std::string> positional;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help") {
            print_help(argv[0]);
            std::exit(0);
        }
        if (argument.rfind("--", 0) == 0) {
            if (index + 1 >= argc) {
                throw std::runtime_error("missing value after " + argument);
            }
            const std::string value = argv[++index];
            if (argument == "--frame-stride") {
                options.frame_stride = parse_positive_integer(value, argument);
            } else if (argument == "--grid-spacing") {
                options.target_grid_spacing =
                    parse_positive_double(value, argument);
            } else if (argument == "--contact-cutoff") {
                options.contact_cutoff = parse_positive_double(value, argument);
            } else if (argument == "--q-max") {
                options.q_max = parse_positive_double(value, argument);
            } else if (argument == "--output") {
                options.metrics_file = value;
            } else if (argument == "--structure-output") {
                options.structure_file = value;
            } else if (argument == "--field-output") {
                options.field_file = value;
            } else if (argument == "--report-output") {
                options.report_file = value;
            } else {
                throw std::runtime_error("unknown option: " + argument);
            }
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() != 2) {
        throw std::runtime_error(
            "expected DATA_OR_TRAJECTORY and INFO_FILE; use --help");
    }
    options.input_file = positional[0];
    options.info_file = positional[1];
    return options;
}

long long parse_integer_line(
    const std::string &line,
    const std::string &description) {
    const std::string value = trim(line);
    std::size_t consumed = 0;
    long long parsed = 0;
    try {
        parsed = std::stoll(value, &consumed);
    } catch (...) {
        throw std::runtime_error("invalid " + description);
    }
    if (consumed != value.size()) {
        throw std::runtime_error("invalid " + description);
    }
    return parsed;
}

std::array<double, 2> parse_bound_line(const std::string &line) {
    std::istringstream input(line);
    std::array<double, 2> bounds{};
    if (!(input >> bounds[0] >> bounds[1])) {
        throw std::runtime_error("invalid BOX BOUNDS line");
    }
    return bounds;
}

std::size_t required_column(
    const std::map<std::string, std::size_t> &columns,
    const std::string &name) {
    const auto found = columns.find(name);
    if (found == columns.end()) {
        throw std::runtime_error("trajectory is missing atom column: " + name);
    }
    return found->second;
}

void parse_numeric_fields(
    const std::string &line,
    std::size_t expected,
    std::array<double, 16> &values) {
    if (expected > values.size()) {
        throw std::runtime_error("trajectory has too many atom columns");
    }
    std::istringstream input(line);
    for (std::size_t index = 0; index < expected; ++index) {
        if (!(input >> values[index])) {
            throw std::runtime_error("incomplete numeric atom line");
        }
    }
    double extra = 0.0;
    if (input >> extra) {
        throw std::runtime_error("too many values in atom line");
    }
}

long long exact_integer(double value, const std::string &field) {
    if (!std::isfinite(value) ||
        std::fabs(value - std::round(value)) > 1.0e-8) {
        throw std::runtime_error("noninteger trajectory " + field);
    }
    return static_cast<long long>(std::llround(value));
}

bool same_box(const Box &left, const Box &right) {
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
        const double scale = std::max(1.0, left.length(dimension));
        if (std::fabs(left.lo[dimension] - right.lo[dimension]) >
                1.0e-10 * scale ||
            std::fabs(left.hi[dimension] - right.hi[dimension]) >
                1.0e-10 * scale ||
            left.periodic[dimension] != right.periodic[dimension]) {
            return false;
        }
    }
    return true;
}

bool read_frame(
    std::ifstream &input,
    const ModelInfo &info,
    long long frame_index,
    Box &reference_box,
    bool &have_reference_box,
    Frame &frame) {
    std::string line;
    while (std::getline(input, line) && trim(line).empty()) {}
    if (!input && trim(line).empty()) return false;
    if (trim(line) != "ITEM: TIMESTEP") {
        throw std::runtime_error(
            "expected ITEM: TIMESTEP before frame " +
            std::to_string(frame_index + 1));
    }
    if (!std::getline(input, line)) {
        throw std::runtime_error("missing trajectory timestep");
    }
    frame = Frame{};
    frame.index = frame_index;
    frame.timestep = parse_integer_line(line, "trajectory timestep");

    if (!std::getline(input, line) || trim(line) != "ITEM: NUMBER OF ATOMS" ||
        !std::getline(input, line)) {
        throw std::runtime_error("invalid NUMBER OF ATOMS block");
    }
    const long long atom_count = parse_integer_line(line, "atom count");
    if (atom_count != info.total_beads) {
        throw std::runtime_error(
            "trajectory atom count does not match info total_beads");
    }

    if (!std::getline(input, line)) {
        throw std::runtime_error("missing BOX BOUNDS header");
    }
    const std::vector<std::string> box_header = split_words(line);
    if (box_header.size() != 6 || box_header[0] != "ITEM:" ||
        box_header[1] != "BOX" || box_header[2] != "BOUNDS") {
        throw std::runtime_error(
            "only orthogonal ITEM: BOX BOUNDS trajectories are supported");
    }
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
        frame.box.periodic[dimension] =
            !box_header[dimension + 3].empty() &&
            box_header[dimension + 3][0] == 'p';
        if (!std::getline(input, line)) {
            throw std::runtime_error("incomplete BOX BOUNDS block");
        }
        const auto bounds = parse_bound_line(line);
        frame.box.lo[dimension] = bounds[0];
        frame.box.hi[dimension] = bounds[1];
        if (frame.box.length(dimension) <= 0.0) {
            throw std::runtime_error("nonpositive trajectory box length");
        }
    }
    if (info.geometry == "bulk" &&
        (!frame.box.periodic[0] || !frame.box.periodic[1] ||
         !frame.box.periodic[2])) {
        throw std::runtime_error("bulk info does not match dump boundaries");
    }
    if (info.geometry == "film" &&
        (!frame.box.periodic[0] || !frame.box.periodic[1] ||
         frame.box.periodic[2])) {
        throw std::runtime_error("film info does not match dump boundaries");
    }
    if (!have_reference_box) {
        reference_box = frame.box;
        have_reference_box = true;
    } else if (!same_box(reference_box, frame.box)) {
        throw std::runtime_error(
            "phase-analysis trajectory box changes between frames");
    }

    if (!std::getline(input, line)) {
        throw std::runtime_error("missing ATOMS header");
    }
    const std::vector<std::string> atom_header = split_words(line);
    if (atom_header.size() < 3 || atom_header[0] != "ITEM:" ||
        atom_header[1] != "ATOMS") {
        throw std::runtime_error("invalid ITEM: ATOMS header");
    }
    const std::size_t field_count = atom_header.size() - 2;
    std::map<std::string, std::size_t> columns;
    for (std::size_t index = 0; index < field_count; ++index) {
        columns[atom_header[index + 2]] = index;
    }
    const std::size_t id_column = required_column(columns, "id");
    const std::size_t molecule_column = required_column(columns, "mol");
    const std::size_t type_column = required_column(columns, "type");
    const std::size_t x_column = required_column(columns, "x");
    const std::size_t y_column = required_column(columns, "y");
    const std::size_t z_column = required_column(columns, "z");

    const auto ends = component_ends(info);
    const long long expected_mps = info.total_mps_repeats;
    const long long expected_dms = info.total_beads - 2 * expected_mps;
    frame.mps.reserve(static_cast<std::size_t>(expected_mps));
    frame.dms.reserve(static_cast<std::size_t>(expected_dms));
    std::vector<unsigned char> seen(
        static_cast<std::size_t>(info.total_beads + 1), 0);
    std::array<double, 16> fields{};

    for (long long atom_line = 0; atom_line < atom_count; ++atom_line) {
        if (!std::getline(input, line)) {
            throw std::runtime_error("trajectory ends inside an atom block");
        }
        parse_numeric_fields(line, field_count, fields);
        const long long atom = exact_integer(fields[id_column], "id");
        const long long molecule =
            exact_integer(fields[molecule_column], "mol");
        const long long type = exact_integer(fields[type_column], "type");
        if (atom < 1 || atom > info.total_beads ||
            seen[static_cast<std::size_t>(atom)] != 0) {
            throw std::runtime_error("invalid or duplicate atom ID in frame");
        }
        seen[static_cast<std::size_t>(atom)] = 1;
        const int component = component_for_molecule(molecule, ends);
        const Vec3 position{
            fields[x_column], fields[y_column], fields[z_column]};
        if (type >= 1 && type <= 3) {
            frame.dms.push_back({atom, molecule, position});
        } else if (type == 4) {
            if (component != kFiller) {
                throw std::runtime_error(
                    "MPS backbone bead lies outside filler molecule range");
            }
            frame.mps.push_back({atom, molecule, position});
        } else if (type == 5) {
            if (component != kFiller) {
                throw std::runtime_error(
                    "MPS pendant bead lies outside filler molecule range");
            }
            ++frame.pendant_count;
        } else {
            throw std::runtime_error(
                "unsupported atom type in phase-analysis trajectory");
        }
    }
    if (static_cast<long long>(frame.mps.size()) != expected_mps ||
        frame.pendant_count != expected_mps ||
        static_cast<long long>(frame.dms.size()) != expected_dms) {
        throw std::runtime_error(
            "DMS/MPS marker counts do not match model-info composition");
    }
    return true;
}

bool input_is_trajectory(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open input file: " + path);
    std::string line;
    while (std::getline(input, line)) {
        if (!trim(line).empty()) return trim(line) == "ITEM: TIMESTEP";
    }
    throw std::runtime_error("input file is empty: " + path);
}

Frame read_data_snapshot(
    const std::string &path,
    const ModelInfo &info) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open LAMMPS data file: " + path);
    }
    Frame frame;
    frame.index = 0;
    frame.timestep = 0;
    frame.box.periodic = {{true, true, info.geometry == "bulk"}};
    long long declared_atoms = -1;
    bool have_x = false;
    bool have_y = false;
    bool have_z = false;
    bool found_atoms = false;
    long long atoms_read = 0;
    std::vector<unsigned char> seen;
    const auto ends = component_ends(info);
    const long long expected_mps = info.total_mps_repeats;
    const long long expected_dms = info.total_beads - 2 * expected_mps;
    frame.mps.reserve(static_cast<std::size_t>(expected_mps));
    frame.dms.reserve(static_cast<std::size_t>(expected_dms));
    std::string line;
    long long line_number = 0;

    while (std::getline(input, line)) {
        ++line_number;
        const std::string clean = trim(line);
        if (!found_atoms) {
            if (clean.empty() || clean.front() == '#') continue;
            std::istringstream header_line(clean);
            long long integer_value = 0;
            std::string label;
            if (header_line >> integer_value >> label) {
                if (label == "atoms") declared_atoms = integer_value;
            }
            std::istringstream bounds_line(clean);
            double lower = 0.0;
            double upper = 0.0;
            std::string lower_label;
            std::string upper_label;
            if (bounds_line >> lower >> upper >> lower_label >> upper_label) {
                if (lower_label == "xlo" && upper_label == "xhi") {
                    frame.box.lo[0] = lower;
                    frame.box.hi[0] = upper;
                    have_x = true;
                } else if (lower_label == "ylo" && upper_label == "yhi") {
                    frame.box.lo[1] = lower;
                    frame.box.hi[1] = upper;
                    have_y = true;
                } else if (lower_label == "zlo" && upper_label == "zhi") {
                    frame.box.lo[2] = lower;
                    frame.box.hi[2] = upper;
                    have_z = true;
                }
            }
            if (clean.rfind("Atoms", 0) == 0) {
                if (declared_atoms < 0 || !have_x || !have_y || !have_z) {
                    throw std::runtime_error(
                        "data-file atom count and box bounds must precede Atoms");
                }
                if (declared_atoms != info.total_beads) {
                    throw std::runtime_error(
                        "data-file atom count does not match info total_beads");
                }
                for (std::size_t dimension = 0; dimension < 3; ++dimension) {
                    if (frame.box.length(dimension) <= 0.0) {
                        throw std::runtime_error(
                            "nonpositive data-file box length");
                    }
                }
                seen.assign(
                    static_cast<std::size_t>(declared_atoms + 1), 0);
                found_atoms = true;
            }
            continue;
        }

        if (atoms_read == declared_atoms) break;
        if (clean.empty() || clean.front() == '#') continue;
        std::istringstream atom_line(clean);
        long long atom = 0;
        long long molecule = 0;
        long long type = 0;
        double charge = 0.0;
        Vec3 position;
        if (!(atom_line >> atom >> molecule >> type >> charge >>
              position.x >> position.y >> position.z)) {
            throw std::runtime_error(
                "invalid atom record at data-file line " +
                std::to_string(line_number));
        }
        (void)charge;
        if (atom < 1 || atom > declared_atoms ||
            seen[static_cast<std::size_t>(atom)] != 0) {
            throw std::runtime_error(
                "invalid or duplicate atom ID in data snapshot");
        }
        seen[static_cast<std::size_t>(atom)] = 1;
        const int component = component_for_molecule(molecule, ends);
        if (type >= 1 && type <= 3) {
            frame.dms.push_back({atom, molecule, position});
        } else if (type == 4) {
            if (component != kFiller) {
                throw std::runtime_error(
                    "MPS backbone bead lies outside filler molecule range");
            }
            frame.mps.push_back({atom, molecule, position});
        } else if (type == 5) {
            if (component != kFiller) {
                throw std::runtime_error(
                    "MPS pendant bead lies outside filler molecule range");
            }
            ++frame.pendant_count;
        } else {
            throw std::runtime_error(
                "unsupported atom type in data snapshot");
        }
        ++atoms_read;
    }
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading LAMMPS data file");
    }
    if (!found_atoms) {
        throw std::runtime_error("missing Atoms section in LAMMPS data file");
    }
    if (atoms_read != declared_atoms) {
        throw std::runtime_error(
            "Atoms section count differs from data-file header");
    }
    if (static_cast<long long>(frame.mps.size()) != expected_mps ||
        frame.pendant_count != expected_mps ||
        static_cast<long long>(frame.dms.size()) != expected_dms) {
        throw std::runtime_error(
            "DMS/MPS marker counts do not match model-info composition");
    }
    return frame;
}

int nearest_power_of_two(int value) {
    int lower = 1;
    while (2 * lower <= value) {
        if (lower > 16384) {
            throw std::runtime_error("composition grid dimension is too large");
        }
        lower *= 2;
    }
    const int upper = 2 * lower;
    return value - lower <= upper - value ? lower : upper;
}

int grid_dimension(double length, double target_spacing) {
    const int requested =
        std::max(4, static_cast<int>(std::ceil(length / target_spacing)));
    const int result = nearest_power_of_two(requested);
    if (result > 256) {
        throw std::runtime_error(
            "grid dimension exceeds 256; increase --grid-spacing");
    }
    return result;
}

int coordinate_index(
    double coordinate,
    double lower,
    double length,
    int cells,
    bool periodic) {
    double reduced = (coordinate - lower) / length;
    if (periodic) {
        reduced -= std::floor(reduced);
    } else {
        reduced = std::max(0.0, std::min(
            std::nextafter(1.0, 0.0), reduced));
    }
    int index = static_cast<int>(std::floor(reduced * cells));
    if (index < 0) index = 0;
    if (index >= cells) index = cells - 1;
    return index;
}

CompositionGrid make_grid(const Frame &frame, double target_spacing) {
    CompositionGrid grid;
    grid.box = frame.box;
    grid.nx = grid_dimension(frame.box.length(0), target_spacing);
    grid.ny = grid_dimension(frame.box.length(1), target_spacing);
    grid.nz = grid_dimension(frame.box.length(2), target_spacing);
    const std::size_t cells =
        static_cast<std::size_t>(grid.nx) *
        static_cast<std::size_t>(grid.ny) *
        static_cast<std::size_t>(grid.nz);
    if (cells > 16777216ULL) {
        throw std::runtime_error(
            "composition grid exceeds 16,777,216 cells; increase --grid-spacing");
    }
    grid.mps.assign(cells, 0);
    grid.dms.assign(cells, 0);
    const auto deposit = [&](const Marker &marker, std::vector<long long> &counts) {
        const int ix = coordinate_index(
            marker.position.x, grid.box.lo[0], grid.box.length(0),
            grid.nx, grid.box.periodic[0]);
        const int iy = coordinate_index(
            marker.position.y, grid.box.lo[1], grid.box.length(1),
            grid.ny, grid.box.periodic[1]);
        const int iz = coordinate_index(
            marker.position.z, grid.box.lo[2], grid.box.length(2),
            grid.nz, grid.box.periodic[2]);
        ++counts[grid.index(ix, iy, iz)];
    };
    for (const Marker &marker : frame.mps) deposit(marker, grid.mps);
    for (const Marker &marker : frame.dms) deposit(marker, grid.dms);
    return grid;
}

SegregationResult calculate_segregation(
    const std::vector<long long> &mps,
    const std::vector<long long> &dms) {
    if (mps.size() != dms.size()) {
        throw std::runtime_error("internal composition-grid size mismatch");
    }
    const long long total_mps =
        std::accumulate(mps.begin(), mps.end(), 0LL);
    const long long total_dms =
        std::accumulate(dms.begin(), dms.end(), 0LL);
    const long long total = total_mps + total_dms;
    if (total_mps <= 0 || total_dms <= 0 || total <= 1) {
        throw std::runtime_error(
            "phase separation requires both MPS and DMS markers");
    }
    const double fraction =
        static_cast<double>(total_mps) / static_cast<double>(total);
    SegregationResult result;
    double observed_weighted_sum = 0.0;
    double random_weighted_sum = 0.0;
    for (std::size_t cell = 0; cell < mps.size(); ++cell) {
        const long long count = mps[cell] + dms[cell];
        if (count <= 0) continue;
        ++result.occupied_cells;
        const double local =
            static_cast<double>(mps[cell]) / static_cast<double>(count);
        const double delta = local - fraction;
        observed_weighted_sum += static_cast<double>(count) * delta * delta;
        random_weighted_sum +=
            fraction * (1.0 - fraction) *
            static_cast<double>(total - count) /
            static_cast<double>(total - 1);
    }
    result.observed_variance =
        observed_weighted_sum / static_cast<double>(total);
    result.random_variance =
        random_weighted_sum / static_cast<double>(total);
    const double fully_separated_variance = fraction * (1.0 - fraction);
    const double denominator =
        fully_separated_variance - result.random_variance;
    result.index = denominator > 0.0
        ? (result.observed_variance - result.random_variance) / denominator
        : 0.0;
    return result;
}

SegregationResult calculate_xy_segregation(const CompositionGrid &grid) {
    const std::size_t cells =
        static_cast<std::size_t>(grid.nx) * static_cast<std::size_t>(grid.ny);
    std::vector<long long> mps(cells, 0);
    std::vector<long long> dms(cells, 0);
    for (int iz = 0; iz < grid.nz; ++iz) {
        for (int iy = 0; iy < grid.ny; ++iy) {
            for (int ix = 0; ix < grid.nx; ++ix) {
                const std::size_t source = grid.index(ix, iy, iz);
                const std::size_t target =
                    static_cast<std::size_t>(iy) *
                        static_cast<std::size_t>(grid.nx) +
                    static_cast<std::size_t>(ix);
                mps[target] += grid.mps[source];
                dms[target] += grid.dms[source];
            }
        }
    }
    return calculate_segregation(mps, dms);
}

void fft_1d(std::vector<std::complex<double>> &values) {
    const std::size_t size = values.size();
    for (std::size_t index = 1, reversed = 0; index < size; ++index) {
        std::size_t bit = size >> 1;
        while (reversed & bit) {
            reversed ^= bit;
            bit >>= 1;
        }
        reversed ^= bit;
        if (index < reversed) std::swap(values[index], values[reversed]);
    }
    for (std::size_t length = 2; length <= size; length <<= 1) {
        const double angle = -2.0 * kPi / static_cast<double>(length);
        const std::complex<double> root(std::cos(angle), std::sin(angle));
        for (std::size_t begin = 0; begin < size; begin += length) {
            std::complex<double> factor(1.0, 0.0);
            for (std::size_t offset = 0; offset < length / 2; ++offset) {
                const std::complex<double> even = values[begin + offset];
                const std::complex<double> odd =
                    values[begin + offset + length / 2] * factor;
                values[begin + offset] = even + odd;
                values[begin + offset + length / 2] = even - odd;
                factor *= root;
            }
        }
    }
}

std::size_t fft_index(int ix, int iy, int iz, int nx, int ny) {
    return (static_cast<std::size_t>(iz) * static_cast<std::size_t>(ny) +
            static_cast<std::size_t>(iy)) *
               static_cast<std::size_t>(nx) +
           static_cast<std::size_t>(ix);
}

void fft_grid(
    std::vector<std::complex<double>> &data,
    int nx,
    int ny,
    int nz) {
    std::vector<std::complex<double>> line;
    line.resize(static_cast<std::size_t>(nx));
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                line[static_cast<std::size_t>(ix)] =
                    data[fft_index(ix, iy, iz, nx, ny)];
            }
            fft_1d(line);
            for (int ix = 0; ix < nx; ++ix) {
                data[fft_index(ix, iy, iz, nx, ny)] =
                    line[static_cast<std::size_t>(ix)];
            }
        }
    }
    line.resize(static_cast<std::size_t>(ny));
    for (int iz = 0; iz < nz; ++iz) {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                line[static_cast<std::size_t>(iy)] =
                    data[fft_index(ix, iy, iz, nx, ny)];
            }
            fft_1d(line);
            for (int iy = 0; iy < ny; ++iy) {
                data[fft_index(ix, iy, iz, nx, ny)] =
                    line[static_cast<std::size_t>(iy)];
            }
        }
    }
    if (nz > 1) {
        line.resize(static_cast<std::size_t>(nz));
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                for (int iz = 0; iz < nz; ++iz) {
                    line[static_cast<std::size_t>(iz)] =
                        data[fft_index(ix, iy, iz, nx, ny)];
                }
                fft_1d(line);
                for (int iz = 0; iz < nz; ++iz) {
                    data[fft_index(ix, iy, iz, nx, ny)] =
                        line[static_cast<std::size_t>(iz)];
                }
            }
        }
    }
}

int signed_mode(int index, int size) {
    return index <= size / 2 ? index : index - size;
}

SpectrumResult calculate_spectrum(
    const CompositionGrid &grid,
    const std::string &geometry,
    double requested_q_max) {
    const bool film = geometry == "film";
    const int nz_fft = film ? 1 : grid.nz;
    const std::size_t cells =
        static_cast<std::size_t>(grid.nx) *
        static_cast<std::size_t>(grid.ny) *
        static_cast<std::size_t>(nz_fft);
    std::vector<long long> mps(cells, 0);
    std::vector<long long> dms(cells, 0);
    if (film) {
        for (int iz = 0; iz < grid.nz; ++iz) {
            for (int iy = 0; iy < grid.ny; ++iy) {
                for (int ix = 0; ix < grid.nx; ++ix) {
                    const std::size_t source = grid.index(ix, iy, iz);
                    const std::size_t target =
                        fft_index(ix, iy, 0, grid.nx, grid.ny);
                    mps[target] += grid.mps[source];
                    dms[target] += grid.dms[source];
                }
            }
        }
    } else {
        mps = grid.mps;
        dms = grid.dms;
    }
    const long long total_mps =
        std::accumulate(mps.begin(), mps.end(), 0LL);
    const long long total_dms =
        std::accumulate(dms.begin(), dms.end(), 0LL);
    const long long total = total_mps + total_dms;
    const double x_mps =
        static_cast<double>(total_mps) / static_cast<double>(total);
    const double x_dms = 1.0 - x_mps;
    if (x_mps <= 0.0 || x_dms <= 0.0) {
        throw std::runtime_error(
            "structure factor requires both MPS and DMS markers");
    }
    std::vector<std::complex<double>> field(cells);
    for (std::size_t cell = 0; cell < cells; ++cell) {
        field[cell] = std::complex<double>(
            x_dms * static_cast<double>(mps[cell]) -
            x_mps * static_cast<double>(dms[cell]), 0.0);
    }
    fft_grid(field, grid.nx, grid.ny, nz_fft);

    const double dx = grid.box.length(0) / static_cast<double>(grid.nx);
    const double dy = grid.box.length(1) / static_cast<double>(grid.ny);
    const double dz = grid.box.length(2) / static_cast<double>(grid.nz);
    double nyquist = std::min(kPi / dx, kPi / dy);
    if (!film) nyquist = std::min(nyquist, kPi / dz);
    SpectrumResult result;
    result.q_max_used = std::isfinite(requested_q_max)
        ? std::min(requested_q_max, 0.95 * nyquist)
        : 0.80 * nyquist;
    double longest_periodic =
        std::max(grid.box.length(0), grid.box.length(1));
    if (!film) longest_periodic =
        std::max(longest_periodic, grid.box.length(2));
    result.dq = 2.0 * kPi / longest_periodic;
    const int bin_count =
        static_cast<int>(std::ceil(result.q_max_used / result.dq)) + 1;
    std::vector<double> q_sum(static_cast<std::size_t>(bin_count), 0.0);
    std::vector<double> s_sum(static_cast<std::size_t>(bin_count), 0.0);
    std::vector<long long> modes(static_cast<std::size_t>(bin_count), 0);
    const double normalization =
        1.0 / (static_cast<double>(total) * x_mps * x_dms);

    for (int iz = 0; iz < nz_fft; ++iz) {
        const int mode_z = film ? 0 : signed_mode(iz, nz_fft);
        const double qz = film ? 0.0 :
            2.0 * kPi * static_cast<double>(mode_z) /
            grid.box.length(2);
        for (int iy = 0; iy < grid.ny; ++iy) {
            const int mode_y = signed_mode(iy, grid.ny);
            const double qy =
                2.0 * kPi * static_cast<double>(mode_y) /
                grid.box.length(1);
            for (int ix = 0; ix < grid.nx; ++ix) {
                const int mode_x = signed_mode(ix, grid.nx);
                if (mode_x == 0 && mode_y == 0 && mode_z == 0) continue;
                const double qx =
                    2.0 * kPi * static_cast<double>(mode_x) /
                    grid.box.length(0);
                const double q = std::sqrt(qx * qx + qy * qy + qz * qz);
                if (q > result.q_max_used) continue;
                int bin = static_cast<int>(std::llround(q / result.dq));
                if (bin <= 0 || bin >= bin_count) continue;
                const double scc =
                    std::norm(field[fft_index(
                        ix, iy, iz, grid.nx, grid.ny)]) *
                    normalization;
                q_sum[static_cast<std::size_t>(bin)] += q;
                s_sum[static_cast<std::size_t>(bin)] += scc;
                ++modes[static_cast<std::size_t>(bin)];
            }
        }
    }
    for (int bin = 1; bin < bin_count; ++bin) {
        const long long count = modes[static_cast<std::size_t>(bin)];
        if (count == 0) continue;
        result.bins.push_back({
            q_sum[static_cast<std::size_t>(bin)] / static_cast<double>(count),
            count,
            s_sum[static_cast<std::size_t>(bin)] / static_cast<double>(count)
        });
    }
    if (result.bins.empty()) {
        throw std::runtime_error(
            "no reciprocal-space modes lie below q_max");
    }
    const std::size_t low_count =
        std::min<std::size_t>(3, result.bins.size());
    for (std::size_t index = 0; index < low_count; ++index) {
        result.low_q_mean += result.bins[index].scc;
    }
    result.low_q_mean /= static_cast<double>(low_count);
    const auto peak = std::max_element(
        result.bins.begin(), result.bins.end(),
        [](const SpectrumBin &left, const SpectrumBin &right) {
            return left.scc < right.scc;
        });
    result.q_peak = peak->q;
    result.peak_scc = peak->scc;
    result.domain_size = result.q_peak > 0.0
        ? 2.0 * kPi / result.q_peak : 0.0;
    result.peak_at_lowest_q = peak == result.bins.begin();
    return result;
}

class DisjointSet {
public:
    explicit DisjointSet(std::size_t size)
        : parent_(size), sizes_(size, 1) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    int find(int value) {
        int root = value;
        while (parent_[static_cast<std::size_t>(root)] != root) {
            root = parent_[static_cast<std::size_t>(root)];
        }
        while (parent_[static_cast<std::size_t>(value)] != value) {
            const int next = parent_[static_cast<std::size_t>(value)];
            parent_[static_cast<std::size_t>(value)] = root;
            value = next;
        }
        return root;
    }

    void unite(int left, int right) {
        left = find(left);
        right = find(right);
        if (left == right) return;
        if (sizes_[static_cast<std::size_t>(left)] <
            sizes_[static_cast<std::size_t>(right)]) {
            std::swap(left, right);
        }
        parent_[static_cast<std::size_t>(right)] = left;
        sizes_[static_cast<std::size_t>(left)] +=
            sizes_[static_cast<std::size_t>(right)];
    }

    int size(int value) {
        return sizes_[static_cast<std::size_t>(find(value))];
    }

private:
    std::vector<int> parent_;
    std::vector<int> sizes_;
};

double minimum_image(double displacement, double length, bool periodic) {
    if (periodic) displacement -= std::round(displacement / length) * length;
    return displacement;
}

ClusterResult calculate_clusters(
    const Frame &frame,
    const ModelInfo &info,
    double cutoff) {
    if (info.oil_chain_count <= 0 || frame.mps.empty()) {
        throw std::runtime_error("MPS cluster analysis requires MPS oil chains");
    }
    const auto ends = component_ends(info);
    const long long first_filler = ends[kCrosslinker] + 1;
    const int chains = static_cast<int>(info.oil_chain_count);
    DisjointSet sets(static_cast<std::size_t>(chains));

    const int nx = std::max(
        1, static_cast<int>(std::floor(frame.box.length(0) / cutoff)));
    const int ny = std::max(
        1, static_cast<int>(std::floor(frame.box.length(1) / cutoff)));
    const int nz = std::max(
        1, static_cast<int>(std::floor(frame.box.length(2) / cutoff)));
    const std::size_t cell_count =
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
        static_cast<std::size_t>(nz);
    std::vector<std::vector<int>> cells(cell_count);
    const auto cell_index = [nx, ny](int ix, int iy, int iz) {
        return (static_cast<std::size_t>(iz) * static_cast<std::size_t>(ny) +
                static_cast<std::size_t>(iy)) *
                   static_cast<std::size_t>(nx) +
               static_cast<std::size_t>(ix);
    };
    std::vector<std::array<int, 3>> marker_cells(frame.mps.size());
    for (std::size_t index = 0; index < frame.mps.size(); ++index) {
        const Marker &marker = frame.mps[index];
        const int ix = coordinate_index(
            marker.position.x, frame.box.lo[0], frame.box.length(0),
            nx, frame.box.periodic[0]);
        const int iy = coordinate_index(
            marker.position.y, frame.box.lo[1], frame.box.length(1),
            ny, frame.box.periodic[1]);
        const int iz = coordinate_index(
            marker.position.z, frame.box.lo[2], frame.box.length(2),
            nz, frame.box.periodic[2]);
        marker_cells[index] = {{ix, iy, iz}};
        cells[cell_index(ix, iy, iz)].push_back(static_cast<int>(index));
    }

    const auto neighbor_coordinate = [](
        int value, int cells_in_dimension, bool periodic, bool &valid) {
        if (periodic) {
            value %= cells_in_dimension;
            if (value < 0) value += cells_in_dimension;
            return value;
        }
        if (value < 0 || value >= cells_in_dimension) {
            valid = false;
            return 0;
        }
        return value;
    };
    std::unordered_set<std::uint64_t> contacts;
    const double cutoff_squared = cutoff * cutoff;
    for (std::size_t index = 0; index < frame.mps.size(); ++index) {
        const auto center = marker_cells[index];
        std::array<std::size_t, 27> neighbor_cells{};
        int unique_cells = 0;
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    bool valid = true;
                    const int ix = neighbor_coordinate(
                        center[0] + dx, nx, frame.box.periodic[0], valid);
                    const int iy = neighbor_coordinate(
                        center[1] + dy, ny, frame.box.periodic[1], valid);
                    const int iz = neighbor_coordinate(
                        center[2] + dz, nz, frame.box.periodic[2], valid);
                    if (!valid) continue;
                    const std::size_t candidate = cell_index(ix, iy, iz);
                    bool duplicate = false;
                    for (int previous = 0; previous < unique_cells; ++previous) {
                        if (neighbor_cells[static_cast<std::size_t>(previous)] ==
                            candidate) {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate) {
                        neighbor_cells[static_cast<std::size_t>(unique_cells++)] =
                            candidate;
                    }
                }
            }
        }
        const Marker &left = frame.mps[index];
        const int left_chain =
            static_cast<int>(left.molecule - first_filler);
        if (left_chain < 0 || left_chain >= chains) {
            throw std::runtime_error(
                "MPS marker molecule lies outside filler chain range");
        }
        for (int neighbor = 0; neighbor < unique_cells; ++neighbor) {
            for (const int candidate_index :
                 cells[neighbor_cells[static_cast<std::size_t>(neighbor)]]) {
                if (candidate_index <= static_cast<int>(index)) continue;
                const Marker &right =
                    frame.mps[static_cast<std::size_t>(candidate_index)];
                if (left.molecule == right.molecule) continue;
                double dx = minimum_image(
                    right.position.x - left.position.x,
                    frame.box.length(0), frame.box.periodic[0]);
                double dy = minimum_image(
                    right.position.y - left.position.y,
                    frame.box.length(1), frame.box.periodic[1]);
                double dz = minimum_image(
                    right.position.z - left.position.z,
                    frame.box.length(2), frame.box.periodic[2]);
                if (dx * dx + dy * dy + dz * dz > cutoff_squared) continue;
                const int right_chain =
                    static_cast<int>(right.molecule - first_filler);
                if (right_chain < 0 || right_chain >= chains) {
                    throw std::runtime_error(
                        "MPS marker molecule lies outside filler chain range");
                }
                sets.unite(left_chain, right_chain);
                const std::uint32_t low =
                    static_cast<std::uint32_t>(
                        std::min(left_chain, right_chain));
                const std::uint32_t high =
                    static_cast<std::uint32_t>(
                        std::max(left_chain, right_chain));
                contacts.insert(
                    (static_cast<std::uint64_t>(low) << 32) |
                    static_cast<std::uint64_t>(high));
            }
        }
    }
    ClusterResult result;
    result.unique_chain_contacts =
        static_cast<long long>(contacts.size());
    for (int chain = 0; chain < chains; ++chain) {
        if (sets.find(chain) == chain) {
            ++result.clusters;
            result.largest_cluster_chains =
                std::max<long long>(
                    result.largest_cluster_chains, sets.size(chain));
        }
    }
    result.largest_cluster_fraction =
        static_cast<double>(result.largest_cluster_chains) /
        static_cast<double>(chains);
    return result;
}

PhaseMetrics analyze_frame(
    const Frame &frame,
    const ModelInfo &info,
    const Options &options,
    long long first_timestep,
    CompositionGrid &grid) {
    grid = make_grid(frame, options.target_grid_spacing);
    PhaseMetrics metrics;
    metrics.frame_index = frame.index;
    metrics.timestep = frame.timestep;
    metrics.time_ns =
        static_cast<double>(frame.timestep - first_timestep) *
        info.timestep_fs / 1.0e6;
    metrics.mps_markers = static_cast<long long>(frame.mps.size());
    metrics.dms_markers = static_cast<long long>(frame.dms.size());
    metrics.mps_fraction =
        static_cast<double>(metrics.mps_markers) /
        static_cast<double>(metrics.mps_markers + metrics.dms_markers);
    metrics.segregation_3d =
        calculate_segregation(grid.mps, grid.dms);
    metrics.segregation_xy = calculate_xy_segregation(grid);
    metrics.spectrum =
        calculate_spectrum(grid, info.geometry, options.q_max);
    metrics.clusters =
        calculate_clusters(frame, info, options.contact_cutoff);
    return metrics;
}

void write_metrics_row(std::ostream &output, const PhaseMetrics &metrics) {
    output << metrics.frame_index << " "
           << metrics.timestep << " "
           << metrics.time_ns << " "
           << metrics.mps_markers << " "
           << metrics.dms_markers << " "
           << metrics.mps_fraction << " "
           << metrics.segregation_3d.index << " "
           << metrics.segregation_3d.observed_variance << " "
           << metrics.segregation_3d.random_variance << " "
           << metrics.segregation_xy.index << " "
           << metrics.spectrum.low_q_mean << " "
           << metrics.spectrum.q_peak << " "
           << metrics.spectrum.domain_size << " "
           << metrics.spectrum.peak_scc << " "
           << (metrics.spectrum.peak_at_lowest_q ? 1 : 0) << " "
           << metrics.clusters.clusters << " "
           << metrics.clusters.largest_cluster_chains << " "
           << metrics.clusters.largest_cluster_fraction << " "
           << metrics.clusters.unique_chain_contacts << "\n";
}

void write_spectrum_rows(
    std::ostream &output,
    const PhaseMetrics &metrics) {
    for (const SpectrumBin &bin : metrics.spectrum.bins) {
        output << metrics.frame_index << " "
               << metrics.timestep << " "
               << metrics.time_ns << " "
               << bin.q << " "
               << bin.modes << " "
               << bin.scc << "\n";
    }
}

void write_final_field(
    const std::string &path,
    const CompositionGrid &grid) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot open final composition field: " + path);
    }
    output << std::scientific << std::setprecision(10);
    const double dx = grid.box.length(0) / static_cast<double>(grid.nx);
    const double dy = grid.box.length(1) / static_cast<double>(grid.ny);
    const double dz = grid.box.length(2) / static_cast<double>(grid.nz);
    for (int iz = 0; iz < grid.nz; ++iz) {
        for (int iy = 0; iy < grid.ny; ++iy) {
            for (int ix = 0; ix < grid.nx; ++ix) {
                const std::size_t cell = grid.index(ix, iy, iz);
                const long long total = grid.mps[cell] + grid.dms[cell];
                const double fraction = total > 0
                    ? static_cast<double>(grid.mps[cell]) /
                        static_cast<double>(total)
                    : 0.0;
                output << ix << " " << iy << " " << iz << " "
                       << grid.box.lo[0] + (static_cast<double>(ix) + 0.5) * dx
                       << " "
                       << grid.box.lo[1] + (static_cast<double>(iy) + 0.5) * dy
                       << " "
                       << grid.box.lo[2] + (static_cast<double>(iz) + 0.5) * dz
                       << " "
                       << grid.mps[cell] << " "
                       << grid.dms[cell] << " "
                       << total << " "
                       << fraction << "\n";
            }
        }
    }
    if (!output) {
        throw std::runtime_error(
            "failed while writing final composition field");
    }
}

void write_report(
    const Options &options,
    const ModelInfo &info,
    bool trajectory_mode,
    long long frames_read,
    long long first_timestep,
    long long last_timestep,
    const CompositionGrid &final_grid,
    const std::vector<PhaseMetrics> &metrics,
    const std::vector<std::string> &warnings) {
    std::ofstream output(options.report_file);
    if (!output) {
        throw std::runtime_error(
            "cannot open phase-analysis report: " + options.report_file);
    }
    const PhaseMetrics &first = metrics.front();
    const PhaseMetrics &last = metrics.back();
    const double dx =
        final_grid.box.length(0) / static_cast<double>(final_grid.nx);
    const double dy =
        final_grid.box.length(1) / static_cast<double>(final_grid.ny);
    const double dz =
        final_grid.box.length(2) / static_cast<double>(final_grid.nz);
    double mean_seg3d = 0.0;
    double mean_segxy = 0.0;
    double mean_low_q = 0.0;
    double mean_largest = 0.0;
    for (const PhaseMetrics &item : metrics) {
        mean_seg3d += item.segregation_3d.index;
        mean_segxy += item.segregation_xy.index;
        mean_low_q += item.spectrum.low_q_mean;
        mean_largest += item.clusters.largest_cluster_fraction;
    }
    const double denominator = static_cast<double>(metrics.size());
    mean_seg3d /= denominator;
    mean_segxy /= denominator;
    mean_low_q /= denominator;
    mean_largest /= denominator;

    output << std::fixed << std::setprecision(8)
           << "DMS/MPS PHASE-SEPARATION ANALYSIS\n"
           << "=================================\n"
           << "Case                        : " << info.case_name << "\n"
           << "Formulation                 : " << info.formulation << "\n"
           << "Geometry                    : " << info.geometry << "\n"
           << "Oil model                   : " << info.oil_model << "\n"
           << "Input mode                  : "
           << (trajectory_mode ? "trajectory" : "final data snapshot") << "\n"
           << "Input file                  : " << options.input_file << "\n"
           << "Model info                  : " << options.info_file << "\n\n"
           << "MARKER DEFINITION\n"
           << "DMS marker                  : one atom of type 1, 2, or 3\n"
           << "MPS marker                  : one type-4 backbone atom per MPS repeat\n"
           << "Type-5 pendant treatment    : validated, not counted a second time\n"
           << "DMS repeats/oil chain       : " << info.dms_repeats_per_chain << "\n"
           << "MPS repeats/oil chain       : " << info.mps_repeats_per_chain << "\n"
           << "Chains with one extra MPS   : " << info.chains_with_extra_mps << "\n"
           << "Total oil MPS repeats       : " << info.total_mps_repeats << "\n"
           << "Oil chains                  : " << info.oil_chain_count << "\n"
           << "Global MPS marker fraction  : " << last.mps_fraction << "\n\n";
    if (trajectory_mode) {
        output << "TRAJECTORY\n"
               << "Frames read                 : " << frames_read
               << " / " << info.expected_frames << " expected from timing\n"
               << "Frames analyzed             : " << metrics.size() << "\n"
               << "Frame stride                : " << options.frame_stride << "\n"
               << "First/last timestep         : " << first_timestep
               << " / " << last_timestep << "\n"
               << "Analyzed time span          : " << last.time_ns << " ns\n\n";
    } else {
        output << "FINAL SNAPSHOT\n"
               << "Snapshots analyzed          : 1\n"
               << "Atoms read                  : " << info.total_beads << "\n"
               << "Time evolution              : unavailable by design\n\n";
    }
    output << "COMPOSITION GRID\n"
           << "Target spacing (A)          : " << options.target_grid_spacing << "\n"
           << "Grid dimensions             : " << final_grid.nx << " "
           << final_grid.ny << " " << final_grid.nz << "\n"
           << "Actual spacing (A)          : " << dx << " " << dy << " " << dz
           << "\n"
           << "3D occupied cells (final)   : "
           << last.segregation_3d.occupied_cells << "\n"
           << "XY occupied cells (final)   : "
           << last.segregation_xy.occupied_cells << "\n\n"
           << "CONCENTRATION STRUCTURE FACTOR\n"
           << "Definition                   : Scc/(x_MPS*x_DMS), ideal random baseline ~1\n"
           << "Dimensionality               : "
           << (info.geometry == "film" ? "2D in-plane q_xy" : "3D radial q")
           << "\n"
           << "q maximum used (1/A)        : " << last.spectrum.q_max_used << "\n"
           << "q shell width (1/A)         : " << last.spectrum.dq << "\n"
           << "Final low-q mean            : " << last.spectrum.low_q_mean << "\n"
           << "Final q peak (1/A)          : " << last.spectrum.q_peak << "\n"
           << "Final 2*pi/q peak (A)       : " << last.spectrum.domain_size << "\n"
           << "Final peak Scc              : " << last.spectrum.peak_scc << "\n"
           << "Peak at lowest q            : "
           << (last.spectrum.peak_at_lowest_q ? "yes" : "no") << "\n\n"
           << "MPS CHAIN CONTACT CLUSTERS\n"
           << "Contact cutoff (A)          : " << options.contact_cutoff << "\n"
           << "Intramolecular contacts     : excluded\n"
           << "Final clusters              : " << last.clusters.clusters << "\n"
           << "Final largest cluster       : "
           << last.clusters.largest_cluster_chains << " chains\n"
           << "Final largest fraction      : "
           << last.clusters.largest_cluster_fraction << "\n"
           << "Final unique chain contacts : "
           << last.clusters.unique_chain_contacts << "\n\n";
    if (trajectory_mode) {
        output << "SUMMARY OVER ANALYZED FRAMES\n"
               << "Mean 3D segregation index   : " << mean_seg3d << "\n"
               << "Mean XY segregation index   : " << mean_segxy << "\n"
               << "Mean low-q Scc              : " << mean_low_q << "\n"
               << "Mean largest-cluster frac.  : " << mean_largest << "\n"
               << "First -> final 3D index     : "
               << first.segregation_3d.index << " -> "
               << last.segregation_3d.index << "\n"
               << "First -> final low-q Scc    : "
               << first.spectrum.low_q_mean << " -> "
               << last.spectrum.low_q_mean << "\n"
               << "First -> final largest frac.: "
               << first.clusters.largest_cluster_fraction << " -> "
               << last.clusters.largest_cluster_fraction << "\n\n";
    } else {
        output << "STATIC SNAPSHOT SUMMARY\n"
               << "3D segregation index        : "
               << last.segregation_3d.index << "\n"
               << "XY segregation index        : "
               << last.segregation_xy.index << "\n"
               << "Low-q Scc                   : "
               << last.spectrum.low_q_mean << "\n"
               << "Largest-cluster fraction    : "
               << last.clusters.largest_cluster_fraction << "\n\n";
    }
    output << "OUTPUTS\n"
           << "Metrics                     : " << options.metrics_file << "\n"
           << "Structure factor            : " << options.structure_file << "\n"
           << "Final composition field     : " << options.field_file << "\n"
           << "Report                      : " << options.report_file << "\n\n"
           << "INTERPRETATION\n"
           << "Positive segregation indices indicate composition variance above the\n"
           << "finite-particle random-label baseline.\n";
    if (trajectory_mode) {
        output << "Growth of low-q Scc, domain spacing, or largest-cluster fraction\n"
               << "indicates aggregation/coarsening.\n";
    } else {
        output << "A final snapshot characterizes morphology but cannot establish\n"
               << "whether the system is equilibrated or still coarsening.\n";
    }
    output
           << "If the peak is at the lowest q, 2*pi/q is box-limited and should be\n"
           << "treated as a lower bound rather than a resolved domain size.\n";
    if (info.geometry == "film") {
        output << "For films, use the XY index and in-plane Scc to distinguish\n"
               << "lateral domains from the z layering measured by z_profile.\n";
    }
    output << "\nCHECKS AND WARNINGS\n";
    if (warnings.empty()) {
        output << "PASS: no consistency warnings.\n";
    } else {
        for (const std::string &warning : warnings) {
            output << "WARNING: " << warning << "\n";
        }
    }
    if (!output) {
        throw std::runtime_error("failed while writing phase-analysis report");
    }
}

}  // namespace

int main(int argc, char **argv) {
    try {
        Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        if (info.total_mps_repeats <= 0 ||
            info.oil_chain_count <= 0) {
            throw std::runtime_error(
                "phase analysis requires a PMPS-containing oil system");
        }
        resolve_output_paths(options, info);
        const bool trajectory_mode = input_is_trajectory(options.input_file);
        std::ofstream metrics_output(options.metrics_file);
        std::ofstream structure_output(options.structure_file);
        if (!metrics_output) {
            throw std::runtime_error(
                "cannot open metrics output: " + options.metrics_file);
        }
        if (!structure_output) {
            throw std::runtime_error(
                "cannot open structure-factor output: " +
                options.structure_file);
        }
        metrics_output << std::scientific << std::setprecision(10);
        structure_output << std::scientific << std::setprecision(10);

        std::vector<std::string> warnings;
        if (trajectory_mode &&
            info.declared_expected_frames != info.expected_frames) {
            warnings.push_back(
                "legacy info declares " +
                std::to_string(info.declared_expected_frames) +
                " frames; timing implies " +
                std::to_string(info.expected_frames) +
                " including timestep zero");
        }
        Box reference_box;
        bool have_reference_box = false;
        Frame pending;
        bool have_pending = false;
        long long frames_read = 0;
        long long first_timestep = 0;
        long long last_timestep = 0;
        long long previous_timestep = -1;
        std::vector<PhaseMetrics> analyzed;
        CompositionGrid final_grid;

        const auto process = [&](const Frame &frame) {
            CompositionGrid grid;
            PhaseMetrics result =
                analyze_frame(frame, info, options, first_timestep, grid);
            write_metrics_row(metrics_output, result);
            write_spectrum_rows(structure_output, result);
            analyzed.push_back(std::move(result));
            final_grid = std::move(grid);
            std::cout << "Analyzed frame " << frame.index
                      << " at timestep " << frame.timestep << "\n";
        };

        if (trajectory_mode) {
            std::ifstream input(options.input_file);
            if (!input) {
                throw std::runtime_error(
                    "cannot open trajectory: " + options.input_file);
            }
            while (true) {
                Frame current;
                if (!read_frame(
                        input, info, frames_read, reference_box,
                        have_reference_box, current)) {
                    break;
                }
                if (frames_read == 0) {
                    first_timestep = current.timestep;
                } else if (current.timestep <= previous_timestep) {
                    throw std::runtime_error(
                        "trajectory timesteps are not increasing");
                }
                previous_timestep = current.timestep;
                if (have_pending &&
                    pending.index % options.frame_stride == 0) {
                    process(pending);
                }
                pending = std::move(current);
                have_pending = true;
                ++frames_read;
            }
            if (!input.good() && !input.eof()) {
                throw std::runtime_error("failed while reading trajectory");
            }
            if (!have_pending) {
                throw std::runtime_error("trajectory contains no frames");
            }
            process(pending);
            last_timestep = pending.timestep;
            if (frames_read != info.expected_frames) {
                warnings.push_back(
                    "trajectory contains " + std::to_string(frames_read) +
                    " frames; timing metadata expects " +
                    std::to_string(info.expected_frames));
            }
        } else {
            Frame snapshot = read_data_snapshot(options.input_file, info);
            frames_read = 1;
            first_timestep = 0;
            last_timestep = 0;
            process(snapshot);
            warnings.push_back(
                "static final-snapshot analysis cannot determine whether "
                "the frozen morphology is equilibrated or still coarsening");
        }
        if (analyzed.back().spectrum.peak_at_lowest_q) {
            warnings.push_back(
                "final Scc peak is at the lowest q; domain size is box-limited");
        }
        write_final_field(options.field_file, final_grid);
        write_report(
            options, info, trajectory_mode, frames_read, first_timestep,
            last_timestep, final_grid, analyzed, warnings);
        if (!metrics_output || !structure_output) {
            throw std::runtime_error("failed while writing phase-analysis tables");
        }
        std::cout << "Wrote phase metrics: " << options.metrics_file << "\n"
                  << "Wrote structure factor: " << options.structure_file << "\n"
                  << "Wrote final composition field: " << options.field_file << "\n"
                  << "Wrote phase report: " << options.report_file << "\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
