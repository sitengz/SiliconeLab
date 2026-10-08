#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "pdms_filler_component.hpp"

namespace {

constexpr double kAvogadroScale = 0.602; // Converts (g/mol)/(g/cm^3) to A^3.
constexpr double kPi = 3.14159265358979323846;
constexpr double kBondLengthAngstrom = pdms_filler::kBondLength;
constexpr double kPlacementSpacing800KAngstrom = 7.5;
constexpr double kFillerZLowerFraction = -0.20;
constexpr double kFillerZUpperFraction = 0.20;
constexpr double kModeratorZLowerFraction = 0.28;
constexpr double kModeratorZUpperFraction = 0.38;
constexpr long long kMsdProductionSteps = 1000000;
constexpr int kMsdDumpEverySteps = 1000;
constexpr double kControlledCuringDensity = 0.5;
constexpr double kControlledCrosslinkProbability = 0.5;
constexpr long long kControlledCompressionStageSteps = 1000000;
constexpr long long kConversionControlledHoldMaximumSteps = 5000000;
constexpr long long kControlledPostCureEquilibrationSteps = 1000000;
constexpr long long kRecommendedTotalBeads = 150000;
constexpr long long kMsdExpectedFrames =
    kMsdProductionSteps / kMsdDumpEverySteps + 1;

struct Settings {
    int n1 = 128, m1 = 900;
    bool strand_length_explicit = false;
    int n2 = 32,  m2 = 0;   // M2 is resolved from stoichiometry.
    int n3 = 5,   m3 = 0;   // Optional five-bead star-like moderators.
    int n4 = 0,   m4 = 0;   // Optional neutral PDMS filler chains.
    int crosslinker_functionality = 4;
    std::string strand_topology = "linear";
    int strand_functionality = 2;
    bool strand_functionality_explicit = false;
    int star_arm_count = 4;
    bool star_arm_count_explicit = false;
    int star_arm_length = 0; // Resolved from N1 when the topology is star.
    int graft_backbone_length = 64;
    int graft_side_chain_length = 16;
    int graft_spacing = 12;
    double graft_functional_fraction = 40.0;
    int graft_side_chain_count = 0;
    int graft_functional_count = 0;
    bool graft_option_explicit = false;
    std::string strand_reactive_distribution = "regular";
    std::uint32_t strand_reactive_seed = 20260810u;
    int stoichiometry_strand_groups = 1;
    int stoichiometry_crosslinker_groups = 1;
    std::string crosslink_distribution = "random";
    std::uint32_t crosslink_seed = 20260722u;
    double crosslink_probability = 0.1;
    double target_conversion_percent = -1.0;
    bool target_conversion_explicit = false;
    double bead_mass = 74.0;
    double density = 0.1;
    double target_density = 0.8;
    // Positive values select a film and specify the nominal 300 K wall-free
    // material thickness. The generated Lz also includes two wall cutoffs.
    double thickness = -1.0;
    std::uint32_t seed = 5489u;
    std::string output = "data.PDMS_elastomer";
    bool output_explicit = false;
    bool filler_length_explicit = false;
    bool filler_weight_explicit = false;
    double filler_weight_percent = -1.0;
    std::uint32_t filler_seed = 20260727u;
    double filler_minimum_separation = 4.5;
    std::string config_file;
};

struct Atom { int id, molecule, type; double charge, x, y, z; };
struct Bond { int id, type, a, b; };
struct Angle { int id, type, a, b, c; };
struct Dihedral { int id, type, a, b, c, d; };
struct Box { double lx, ly, lz; };

struct System {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
};

int parse_int(const std::string& value, const std::string& option) {
    try { size_t used = 0; int result = std::stoi(value, &used); if (used != value.size()) throw std::invalid_argument(""); return result; }
    catch (...) { throw std::runtime_error("Invalid integer for " + option + ": " + value); }
}

double parse_double(const std::string& value, const std::string& option) {
    try { size_t used = 0; double result = std::stod(value, &used); if (used != value.size()) throw std::invalid_argument(""); return result; }
    catch (...) { throw std::runtime_error("Invalid number for " + option + ": " + value); }
}

void parse_stoichiometry(const std::string& value, Settings& settings) {
    const std::size_t separator = value.find(':');
    if (separator == std::string::npos || value.find(':', separator + 1) != std::string::npos)
        throw std::runtime_error(
            "--stoichiometry must use STRAND:CROSSLINKER form, for example 1:1");
    settings.stoichiometry_strand_groups =
        parse_int(value.substr(0, separator), "--stoichiometry");
    settings.stoichiometry_crosslinker_groups =
        parse_int(value.substr(separator + 1), "--stoichiometry");
    if (settings.stoichiometry_strand_groups <= 0 ||
        settings.stoichiometry_crosslinker_groups <= 0)
        throw std::runtime_error("--stoichiometry values must be positive integers");
}

std::string trim(std::string value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string option_name(std::string key) {
    key = trim(key);
    if (key.rfind("--", 0) == 0) key.erase(0, 2);
    std::replace(key.begin(), key.end(), '_', '-');
    return "--" + key;
}

void apply_option(Settings& s, const std::string& option,
                  const std::string& value) {
    if      (option == "--strand-length") {
        s.n1 = parse_int(value, option);
        s.strand_length_explicit = true;
    }
    else if (option == "--strand-count") s.m1 = parse_int(value, option);
    else if (option == "--strand-topology") s.strand_topology = value;
    else if (option == "--strand-functionality") {
        s.strand_functionality = parse_int(value, option);
        s.strand_functionality_explicit = true;
    }
    else if (option == "--strand-arm-count") {
        s.star_arm_count = parse_int(value, option);
        s.star_arm_count_explicit = true;
    }
    else if (option == "--backbone-length") {
        s.graft_backbone_length = parse_int(value, option);
        s.graft_option_explicit = true;
    }
    else if (option == "--side-chain-length") {
        s.graft_side_chain_length = parse_int(value, option);
        s.graft_option_explicit = true;
    }
    else if (option == "--graft-spacing") {
        s.graft_spacing = parse_int(value, option);
        s.graft_option_explicit = true;
    }
    else if (option == "--graft-functional-fraction") {
        s.graft_functional_fraction = parse_double(value, option);
        s.graft_option_explicit = true;
    }
    else if (option == "--strand-reactive-distribution")
        s.strand_reactive_distribution = value;
    else if (option == "--strand-reactive-seed")
        s.strand_reactive_seed =
            static_cast<std::uint32_t>(parse_int(value, option));
    else if (option == "--crosslinker-length") s.n2 = parse_int(value, option);
    else if (option == "--stoichiometry") parse_stoichiometry(value, s);
    else if (option == "--moderator-count") s.m3 = parse_int(value, option);
    else if (option == "--filler-length") {
        s.n4 = parse_int(value, option);
        s.filler_length_explicit = true;
    }
    else if (option == "--filler-wt") {
        s.filler_weight_percent = parse_double(value, option);
        s.filler_weight_explicit = true;
    }
    else if (option == "--filler-seed")
        s.filler_seed = static_cast<std::uint32_t>(parse_int(value, option));
    else if (option == "--filler-min-separation")
        s.filler_minimum_separation = parse_double(value, option);
    else if (option == "--functionality")
        s.crosslinker_functionality = parse_int(value, option);
    else if (option == "--crosslink-distribution")
        s.crosslink_distribution = value;
    else if (option == "--crosslink-seed")
        s.crosslink_seed = static_cast<std::uint32_t>(parse_int(value, option));
    else if (option == "--target-conversion") {
        s.target_conversion_percent = parse_double(value, option);
        s.target_conversion_explicit = true;
    }
    else if (option == "--mass") s.bead_mass = parse_double(value, option);
    else if (option == "--density") s.density = parse_double(value, option);
    else if (option == "--target-density")
        s.target_density = parse_double(value, option);
    else if (option == "--thickness") s.thickness = parse_double(value, option);
    else if (option == "--seed")
        s.seed = static_cast<std::uint32_t>(parse_int(value, option));
    else if (option == "--output") {
        s.output = value;
        s.output_explicit = true;
    }
    else throw std::runtime_error("Unknown option: " + option);
}

void read_config_file(Settings& settings, const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open config file: " + path);
    settings.config_file = path;

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
            apply_option(settings, key, value);
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Config line " + std::to_string(line_number) + ": " +
                error.what());
        }
    }
    if (!input.good() && !input.eof())
        throw std::runtime_error("Failed while reading config file: " + path);
}

void print_help(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Generic coarse-grained PDMS elastomer generator.\n\n"
        << "Components are written in this molecule-ID order:\n"
        << "  1. linear, ring, star, or grafted-backbone strands\n"
        << "  2. functional short-chain crosslinkers\n"
        << "  3. optional five-bead star moderators\n"
        << "  4. optional neutral PDMS filler chains\n\n"
        << "  --strand-length N       beads per linear/ring strand, or per star arm\n"
        << "                          (default: 128; star-arm default: 32)\n"
        << "  --strand-count M        number of strands (default: 900)\n"
        << "  --strand-topology TYPE  linear, ring, star, or grafted (default: linear)\n"
        << "  --strand-functionality F reactive sites; derived from arms for stars\n"
        << "  --strand-arm-count A    star arms: 3, 4, 6, or 8 (default: 4)\n"
        << "  --backbone-length N     grafted backbone beads (default: 64)\n"
        << "  --side-chain-length N   beads per grafted side chain (default: 16)\n"
        << "  --graft-spacing N       ungrafted beads between grafts; 0 grafts every bead\n"
        << "                          (default: 12)\n"
        << "  --graft-functional-fraction X\n"
        << "                          percent of side-chain ends that are functional\n"
        << "                          (default: 40)\n"
        << "  --strand-reactive-distribution MODE\n"
        << "                          ring sites: regular or random (default: regular)\n"
        << "  --strand-reactive-seed N random ring-site seed (default: 20260810)\n"
        << "  --crosslinker-length N  beads per crosslinker (default: 32)\n"
        << "  --functionality F       reactive sites per crosslinker (default: 4)\n"
        << "  --stoichiometry A:B     strand-group : crosslinker-group ratio\n"
        << "                          (default: 1:1; gives strand:crosslinker M = 2:1)\n"
        << "  --moderator-count M     optional five-bead star moderators (default: 0)\n"
        << "  --config FILE           read key = value settings; CLI values override file\n"
        << "\nOptional component-4 PDMS filler:\n"
        << "  --filler-length N repeat units per filler chain\n"
        << "  --filler-wt X     filler weight percent of the complete model\n"
        << "  --filler-seed N   filler conformation/packing seed (default: 20260727)\n"
        << "  --filler-min-separation X  minimum filler/moderator-to-other distance in A"
           " (default: 4.5)\n"
        << "\nNetwork and box controls:\n"
        << "  --crosslink-distribution MODE\n"
        << "                     reactive-site placement: regular or random (default: random)\n"
        << "  --crosslink-seed N random-site seed (default: 20260722)\n"
        << "  --target-conversion X percent\n"
        << "                     stop bond/create at X percent of the stoichiometric\n"
        << "                     maximum; omit to retain the time-controlled workflow\n"
        << "  --mass X          bead mass in g/mol (default: 74)\n"
        << "  --density X       initial packing density in g/cm^3 (default: 0.1)\n"
        << "  --target-density X density after scripted compression (default: 0.8)\n"
        << "  --thickness X     nominal 300 K wall-free film thickness in angstrom;\n"
        << "                    generated Lz = X + 2*wall cutoff; omit for bulk\n"
        << "  --seed N          placement and initial-velocity seed (default: 5489)\n"
        << "  --output FILE     override the automatically generated data filename\n"
        << "                    a case folder is created beside this path\n"
        << "  --help             show this help\n";
}

int star_center_count(int arm_count) {
    if (arm_count == 3 || arm_count == 4) return 1;
    if (arm_count == 6) return 2;
    if (arm_count == 8) return 3;
    throw std::runtime_error("Star arm count must be 3, 4, 6, or 8");
}

void resolve_strand_architecture(Settings& s) {
    if (s.strand_topology == "grafted") {
        if (s.star_arm_count_explicit)
            throw std::runtime_error(
                "--strand-arm-count requires --strand-topology star");
        if (s.strand_length_explicit)
            throw std::runtime_error(
                "Use --backbone-length and --side-chain-length for grafted strands");
        if (s.graft_backbone_length <= 0 || s.graft_side_chain_length <= 0)
            throw std::runtime_error(
                "Grafted backbone and side-chain lengths must be positive");
        if (s.graft_spacing < 0)
            throw std::runtime_error("Graft spacing cannot be negative");
        if (s.graft_functional_fraction <= 0.0 ||
            s.graft_functional_fraction > 100.0)
            throw std::runtime_error(
                "Graft functional fraction must be greater than 0 and at most 100 percent");
        const long long interval = 1LL + s.graft_spacing;
        s.graft_side_chain_count = static_cast<int>(
            (s.graft_backbone_length + interval - 1) / interval);
        s.graft_functional_count = std::max(
            1, static_cast<int>(std::lround(
                s.graft_side_chain_count *
                s.graft_functional_fraction / 100.0)));
        s.graft_functional_count = std::min(
            s.graft_functional_count, s.graft_side_chain_count);
        if (s.strand_functionality_explicit &&
            s.strand_functionality != s.graft_functional_count)
            throw std::runtime_error(
                "Grafted strand functionality is derived from the functional fraction");
        const long long total_beads = s.graft_backbone_length +
            1LL*s.graft_side_chain_count*s.graft_side_chain_length;
        if (total_beads > 100000000LL)
            throw std::runtime_error("Grafted strand bead count is too large");
        s.n1 = static_cast<int>(total_beads);
        s.strand_functionality = s.graft_functional_count;
        return;
    }

    if (s.strand_topology != "star") {
        if (s.star_arm_count_explicit)
            throw std::runtime_error(
                "--strand-arm-count requires --strand-topology star");
        if (s.graft_option_explicit)
            throw std::runtime_error(
                "Grafted-backbone options require --strand-topology grafted");
        return;
    }

    if (s.graft_option_explicit)
        throw std::runtime_error(
            "Grafted-backbone options require --strand-topology grafted");

    const int centers = star_center_count(s.star_arm_count);
    if (!s.strand_length_explicit) s.n1 = 32;
    if (s.n1 <= 0)
        throw std::runtime_error("Star arm length must be positive");
    if (s.strand_functionality_explicit &&
        s.strand_functionality != s.star_arm_count)
        throw std::runtime_error(
            "Star strand functionality must equal its arm count");
    const long long total_beads =
        1LL * s.star_arm_count * s.n1 + centers;
    if (total_beads > 100000000LL)
        throw std::runtime_error("Star strand bead count is too large");
    s.star_arm_length = s.n1;
    s.n1 = static_cast<int>(total_beads);
    s.strand_functionality = s.star_arm_count;
}

Settings parse_args(int argc, char** argv) {
    Settings s;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--help") {
            print_help(argv[0]);
            std::exit(0);
        }
    }

    bool have_config = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) != "--config") continue;
        if (have_config)
            throw std::runtime_error("--config may be supplied only once");
        if (i + 1 >= argc)
            throw std::runtime_error("Missing value after --config");
        read_config_file(s, argv[++i]);
        have_config = true;
    }

    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (i + 1 >= argc) throw std::runtime_error("Missing value after " + option);
        const std::string value = argv[++i];
        if (option != "--config") apply_option(s, option, value);
    }
    return s;
}

std::string filename_number(double value) {
    std::ostringstream text;
    text << std::fixed << std::setprecision(4) << value;
    std::string result = text.str();
    while (!result.empty() && result.back() == '0') result.pop_back();
    if (!result.empty() && result.back() == '.') result.pop_back();
    std::replace(result.begin(), result.end(), '.', 'p');
    return result;
}

double filler_chain_mass(const Settings& s) {
    return pdms_filler::chain_mass(s.n4, s.bead_mass);
}

long long filler_beads(const Settings& s) {
    return 1LL * s.m4 * s.n4;
}

void resolve_filler_composition(Settings& s) {
    const bool requested = s.filler_length_explicit || s.filler_weight_explicit;
    if (!requested) {
        s.n4 = 0;
        s.m4 = 0;
        return;
    }
    if (!s.filler_length_explicit || !s.filler_weight_explicit)
        throw std::runtime_error(
            "PDMS filler requires both --filler-length and --filler-wt");
    if (s.n4 <= 0)
        throw std::runtime_error("--filler-length must be positive");
    if (s.filler_weight_percent <= 0.0 || s.filler_weight_percent >= 100.0)
        throw std::runtime_error(
            "--filler-wt must be greater than 0 and less than 100");
    if (s.filler_minimum_separation <= 0.0 ||
        s.filler_minimum_separation >= 15.0)
        throw std::runtime_error(
            "--filler-min-separation must be greater than 0 and less than 15 A");

    const long long base_beads =
        1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n3*s.m3;
    const double base_mass = base_beads * s.bead_mass;
    const double fraction = s.filler_weight_percent / 100.0;
    const double desired_filler_mass = base_mass * fraction / (1.0 - fraction);
    const double requested_chains = desired_filler_mass / filler_chain_mass(s);
    if (!std::isfinite(requested_chains) ||
        requested_chains > std::numeric_limits<int>::max())
        throw std::runtime_error(
            "Requested filler chain count is too large to represent");
    s.m4 = std::max(
        1, static_cast<int>(std::lround(requested_chains)));

    if (!s.output_explicit) {
        s.output = std::string("data.PDMS_elastomer_filler_N") +
            std::to_string(s.n4) + "_" +
            filename_number(s.filler_weight_percent) + "wt";
    }
}

void apply_crosslinker_stoichiometry(Settings& s) {
    if (s.crosslinker_functionality < 3)
        throw std::runtime_error("Cross-linker functionality must be at least 3 to form a network");
    const long long strand_side_groups =
        1LL * s.strand_functionality * s.m1;
    const long long scaled_crosslinker_groups =
        strand_side_groups * s.stoichiometry_crosslinker_groups;
    const long long denominator =
        1LL * s.crosslinker_functionality * s.stoichiometry_strand_groups;
    if (scaled_crosslinker_groups % denominator != 0) {
        std::ostringstream message;
        message << "Exact stoichiometry is impossible: " << strand_side_groups
                << " strand functional groups at " << s.stoichiometry_strand_groups
                << ':' << s.stoichiometry_crosslinker_groups
                << " cannot be represented by whole crosslinkers with functionality "
                << s.crosslinker_functionality;
        throw std::runtime_error(message.str());
    }
    const long long molecule_count = scaled_crosslinker_groups / denominator;
    if (molecule_count > 100000000LL)
        throw std::runtime_error("Stoichiometric M2 is too large");
    s.m2 = static_cast<int>(molecule_count);
}

long long strand_functional_groups(const Settings& s) {
    return 1LL * s.strand_functionality * s.m1;
}

long long crosslinker_functional_groups(const Settings& s) {
    return 1LL * s.crosslinker_functionality * s.m2;
}

long long stoichiometric_maximum_bonds(const Settings& s) {
    return std::min(strand_functional_groups(s),
                    crosslinker_functional_groups(s));
}

bool conversion_control_enabled(const Settings& s) {
    return s.target_conversion_explicit;
}

long long target_new_bonds(const Settings& s) {
    if (!conversion_control_enabled(s)) return -1;
    return static_cast<long long>(std::floor(
        stoichiometric_maximum_bonds(s) *
        s.target_conversion_percent / 100.0 + 1.0e-12));
}

void apply_geometry_filename(Settings& s) {
    if (s.thickness > 0.0 && !s.output_explicit)
        s.output += "_film_H" + filename_number(s.thickness);
}

void report_composition(const Settings& s) {
    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    const long long component_beads[] = {
        1LL*s.n1*s.m1,
        1LL*s.n2*s.m2,
        1LL*s.n3*s.m3,
        filler_beads(s)
    };
    const double component_masses[] = {
        component_beads[0] * s.bead_mass,
        component_beads[1] * s.bead_mass,
        component_beads[2] * s.bead_mass,
        s.m4 * filler_chain_mass(s)
    };
    long long total_beads = 0;
    long long total_molecules = 0;
    double total_mass = 0.0;
    for (size_t i = 0; i < 4; ++i) {
        total_beads += component_beads[i];
        total_mass += component_masses[i];
    }
    for (int count : m) total_molecules += count;
    std::cerr << "Composition:\n";
    for (size_t i = 0; i < 4; ++i) {
        const double realized = total_mass == 0.0
            ? 0.0 : 100.0 * component_masses[i] / total_mass;
        const double mole_fraction = total_molecules == 0 ? 0.0 : 100.0 * m[i] / total_molecules;
        std::cerr << "  component " << i + 1 << ": N=" << n[i] << ", M=" << m[i]
                  << ", beads=" << component_beads[i]
                  << ", mole%=" << std::fixed << std::setprecision(4) << mole_fraction
                  << ", realized wt%=" << realized;
        if (i == 3 && s.m4 > 0)
            std::cerr << ", filler=PDMS"
                      << ", requested filler wt%=" << s.filler_weight_percent;
        std::cerr << '\n';
    }
    std::cerr << "  total beads=" << total_beads
              << ", total mass=" << total_mass << " g/mol-equivalent\n";
    if (total_beads > kRecommendedTotalBeads)
        std::cerr << "  warning: total beads exceed the recommended "
                  << kRecommendedTotalBeads
                  << "-bead balance between simulation cost and size effects; "
                     "generation will continue\n";
    if (s.m3 > 0)
        std::cerr << "  moderator functional groups=" << 4LL*s.m3
                  << " (extra; excluded from stoichiometry)\n";
    if (conversion_control_enabled(s))
        std::cerr << "  target conversion=" << s.target_conversion_percent
                  << "% of " << stoichiometric_maximum_bonds(s)
                  << " stoichiometric bonds; target new bonds="
                  << target_new_bonds(s) << '\n';
}

void validate(const Settings& s) {
    const int values[] = {s.n1, s.n2, s.n3, s.n4, s.m1, s.m2, s.m3, s.m4};
    if (std::any_of(std::begin(values), std::end(values), [](int x) { return x < 0; }))
        throw std::runtime_error("N and M values cannot be negative");
    if (s.n1 <= 0 || s.m1 <= 0)
        throw std::runtime_error("Strand length and strand count must be positive");
    if (s.strand_topology != "linear" && s.strand_topology != "ring" &&
        s.strand_topology != "star" && s.strand_topology != "grafted")
        throw std::runtime_error(
            "--strand-topology must be linear, ring, star, or grafted");
    const int minimum_strand_functionality =
        s.strand_topology == "grafted" ? 1 : 2;
    if (s.strand_functionality < minimum_strand_functionality ||
        s.strand_functionality > s.n1)
        throw std::runtime_error(
            "Strand functionality is outside the valid range");
    if (s.strand_reactive_distribution != "regular" &&
        s.strand_reactive_distribution != "random")
        throw std::runtime_error(
            "--strand-reactive-distribution must be regular or random");
    if (s.strand_topology == "linear" && s.strand_functionality != 2)
        throw std::runtime_error(
            "Linear strands are bifunctional; use --strand-functionality 2");
    if (s.strand_topology != "ring" &&
        s.strand_reactive_distribution != "regular")
        throw std::runtime_error(
            "Random strand reactive sites require --strand-topology ring");
    if (s.strand_topology == "ring" && s.n1 < 4)
        throw std::runtime_error("Ring strands require at least 4 beads");
    if (s.strand_topology == "ring" &&
        s.strand_reactive_distribution == "regular" &&
        s.n1 % s.strand_functionality != 0)
        throw std::runtime_error(
            "Regular ring reactive sites require strand length divisible by functionality");
    if (s.strand_topology == "star") {
        if (s.star_arm_length <= 0)
            throw std::runtime_error("Star arm length must be positive");
        if (s.strand_functionality != s.star_arm_count)
            throw std::logic_error(
                "Internal error: star functionality does not match arm count");
        star_center_count(s.star_arm_count);
    }
    if (s.strand_topology == "grafted") {
        if (s.graft_backbone_length <= 0 || s.graft_side_chain_length <= 0 ||
            s.graft_side_chain_count <= 0 || s.graft_functional_count <= 0)
            throw std::logic_error(
                "Internal error: invalid resolved grafted architecture");
        if (s.strand_functionality != s.graft_functional_count)
            throw std::logic_error(
                "Internal error: grafted functionality mismatch");
    }
    if (s.n3 != 5 && s.m3 != 0)
        throw std::runtime_error("The implemented moderator is a five-bead star");
    if (s.crosslinker_functionality < 3 || s.crosslinker_functionality > 16 || s.crosslinker_functionality > s.n2)
        throw std::runtime_error("Cross-linker functionality must be between 3 and min(16, N2)");
    if (s.crosslink_distribution != "regular" && s.crosslink_distribution != "random")
        throw std::runtime_error("--crosslink-distribution must be regular or random");
    if (s.crosslink_distribution == "regular" && s.crosslinker_functionality > (s.n2 + 1) / 2)
        throw std::runtime_error("Regular placement at 1,3,5,... requires functionality <= ceil(N2/2)");
    if (s.bead_mass <= 0 || s.density <= 0 || s.target_density <= 0)
        throw std::runtime_error("Mass and densities must be positive");
    if (conversion_control_enabled(s) &&
        (!std::isfinite(s.target_conversion_percent) ||
         s.target_conversion_percent <= 0.0 ||
         s.target_conversion_percent > 100.0))
        throw std::runtime_error(
            "--target-conversion must be greater than 0 and at most 100 percent");
    if (conversion_control_enabled(s) && target_new_bonds(s) < 1)
        throw std::runtime_error(
            "--target-conversion produces fewer than one target bond for this composition");
    if (conversion_control_enabled(s) &&
        (s.density > kControlledCuringDensity ||
         s.target_density < kControlledCuringDensity))
        throw std::runtime_error(
            "Conversion-controlled curing requires initial density <= 0.5 "
            "and target density >= 0.5 g/cm^3");
    if (s.thickness == 0.0)
        throw std::runtime_error("--thickness must be positive; omit it for the cubic bulk system");
    const long long total_beads =
        1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n3*s.m3 + filler_beads(s);
    if (total_beads <= 0) throw std::runtime_error("The model must contain at least one bead");
    if (s.output.empty()) throw std::runtime_error("Output filename cannot be empty");
}

std::set<int> regular_sites(int functionality) {
    std::set<int> sites;
    for (int i = 0; i < functionality; ++i) sites.insert(1 + 2*i);
    return sites;
}

std::set<int> random_sites(int bead_count, int functionality, std::uint32_t seed) {
    std::vector<int> candidates(static_cast<size_t>(bead_count));
    for (int i = 0; i < bead_count; ++i) candidates[static_cast<size_t>(i)] = i + 1;
    std::mt19937 rng(seed);
    std::shuffle(candidates.begin(), candidates.end(), rng);
    return std::set<int>(candidates.begin(), candidates.begin() + functionality);
}

std::vector<int> graft_backbone_sites(const Settings& s) {
    std::vector<int> sites;
    const int interval = s.graft_spacing + 1;
    for (int site = 0; site < s.graft_backbone_length; site += interval)
        sites.push_back(site);
    return sites;
}

std::set<int> graft_functional_side_indices(const Settings& s) {
    std::set<int> indices;
    for (int functional = 0; functional < s.graft_functional_count;
         ++functional) {
        int index = static_cast<int>(std::floor(
            (functional + 0.5) * s.graft_side_chain_count /
            s.graft_functional_count));
        index = std::min(index, s.graft_side_chain_count - 1);
        indices.insert(index);
    }
    if (static_cast<int>(indices.size()) != s.graft_functional_count)
        throw std::logic_error(
            "Internal error: grafted functional-site selection is not unique");
    return indices;
}

std::set<int> strand_sites(const Settings& s) {
    if (s.strand_topology == "linear") return {1, s.n1};
    if (s.strand_topology == "grafted") {
        std::set<int> sites;
        for (int side : graft_functional_side_indices(s))
            sites.insert(
                s.graft_backbone_length +
                (side + 1)*s.graft_side_chain_length);
        return sites;
    }
    if (s.strand_topology == "star") {
        std::set<int> sites;
        const int centers = star_center_count(s.star_arm_count);
        for (int arm = 1; arm <= s.star_arm_count; ++arm)
            sites.insert(centers + arm * s.star_arm_length);
        return sites;
    }
    if (s.strand_reactive_distribution == "random")
        return random_sites(s.n1, s.strand_functionality,
                            s.strand_reactive_seed);
    std::set<int> sites;
    const int interval = s.n1 / s.strand_functionality;
    for (int i = 0; i < s.strand_functionality; ++i)
        sites.insert(1 + i * interval);
    return sites;
}

std::string strand_reactive_distribution_name(const Settings& s) {
    if (s.strand_topology == "ring")
        return s.strand_reactive_distribution;
    if (s.strand_topology == "star") return "arm_ends";
    if (s.strand_topology == "grafted")
        return "selected_side_chain_ends";
    return "ends";
}

std::string grafted_architecture_label(const Settings& s) {
    if (s.strand_topology != "grafted") return "not_applicable";
    return s.graft_spacing == 0 ? "dense_bottlebrush" : "comb_like";
}

std::set<int> crosslinker_sites(const Settings& s) {
    if (s.crosslink_distribution == "regular")
        return regular_sites(s.crosslinker_functionality);
    return random_sites(s.n2, s.crosslinker_functionality, s.crosslink_seed);
}

void add_linear_topology(System& sys, int first_atom, int bead_count) {
    for (int i = 0; i + 1 < bead_count; ++i)
        sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 1, first_atom + i, first_atom + i + 1});
    for (int i = 0; i + 2 < bead_count; ++i)
        sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1, first_atom + i, first_atom + i + 1, first_atom + i + 2});
    for (int i = 0; i + 3 < bead_count; ++i)
        sys.dihedrals.push_back({static_cast<int>(sys.dihedrals.size()) + 1, 1, first_atom + i, first_atom + i + 1, first_atom + i + 2, first_atom + i + 3});
}

void add_ring_topology(System& sys, int first_atom, int bead_count) {
    const auto atom = [first_atom, bead_count](int offset) {
        return first_atom + (offset % bead_count);
    };
    for (int i = 0; i < bead_count; ++i)
        sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 1,
                             atom(i), atom(i + 1)});
    for (int i = 0; i < bead_count; ++i)
        sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1,
                              atom(i), atom(i + 1), atom(i + 2)});
    for (int i = 0; i < bead_count; ++i)
        sys.dihedrals.push_back({static_cast<int>(sys.dihedrals.size()) + 1, 1,
                                 atom(i), atom(i + 1), atom(i + 2), atom(i + 3)});
}

struct StarArchitecture {
    int center_count = 0;
    std::vector<int> attachment_centers;
    std::vector<std::pair<int, int>> edges;
};

struct GraftedArchitecture {
    std::vector<int> graft_sites;
    std::vector<std::pair<int, int>> edges;
};

GraftedArchitecture grafted_architecture(const Settings& s) {
    GraftedArchitecture architecture;
    architecture.graft_sites = graft_backbone_sites(s);
    if (static_cast<int>(architecture.graft_sites.size()) !=
        s.graft_side_chain_count)
        throw std::logic_error(
            "Internal error: resolved graft count does not match graft sites");
    for (int bead = 0; bead + 1 < s.graft_backbone_length; ++bead)
        architecture.edges.push_back({bead, bead + 1});
    for (int side = 0; side < s.graft_side_chain_count; ++side) {
        const int first = s.graft_backbone_length +
                          side*s.graft_side_chain_length;
        architecture.edges.push_back({
            architecture.graft_sites[static_cast<std::size_t>(side)], first
        });
        for (int bead = 1; bead < s.graft_side_chain_length; ++bead)
            architecture.edges.push_back({first + bead - 1, first + bead});
    }
    return architecture;
}

StarArchitecture star_architecture(const Settings& s) {
    StarArchitecture architecture;
    architecture.center_count = star_center_count(s.star_arm_count);
    architecture.attachment_centers.reserve(
        static_cast<std::size_t>(s.star_arm_count));
    if (s.star_arm_count == 3 || s.star_arm_count == 4) {
        architecture.attachment_centers.assign(
            static_cast<std::size_t>(s.star_arm_count), 0);
    } else if (s.star_arm_count == 6) {
        architecture.attachment_centers = {0, 0, 0, 1, 1, 1};
    } else {
        architecture.attachment_centers = {0, 0, 0, 1, 1, 2, 2, 2};
    }

    for (int center = 0; center + 1 < architecture.center_count; ++center)
        architecture.edges.push_back({center, center + 1});
    for (int arm = 0; arm < s.star_arm_count; ++arm) {
        const int first = architecture.center_count + arm * s.star_arm_length;
        architecture.edges.push_back({
            architecture.attachment_centers[static_cast<std::size_t>(arm)],
            first
        });
        for (int bead = 1; bead < s.star_arm_length; ++bead)
            architecture.edges.push_back({first + bead - 1, first + bead});
    }
    return architecture;
}

std::vector<std::vector<int>> graph_adjacency(
    int bead_count,
    const std::vector<std::pair<int, int>>& edges
) {
    std::vector<std::vector<int>> adjacency(
        static_cast<std::size_t>(bead_count));
    for (const auto& edge : edges) {
        adjacency[static_cast<std::size_t>(edge.first)].push_back(edge.second);
        adjacency[static_cast<std::size_t>(edge.second)].push_back(edge.first);
    }
    return adjacency;
}

struct TopologyCounts {
    long long bonds = 0;
    long long angles = 0;
    long long dihedrals = 0;
};

TopologyCounts graph_topology_counts(
    int bead_count,
    const std::vector<std::pair<int, int>>& edges
) {
    const std::vector<std::vector<int>> adjacency =
        graph_adjacency(bead_count, edges);
    TopologyCounts counts;
    counts.bonds = static_cast<long long>(edges.size());
    for (const auto& neighbors : adjacency) {
        const long long degree = static_cast<long long>(neighbors.size());
        counts.angles += degree * (degree - 1) / 2;
    }
    for (const auto& edge : edges) {
        counts.dihedrals +=
            (static_cast<long long>(adjacency[static_cast<std::size_t>(edge.first)].size()) - 1) *
            (static_cast<long long>(adjacency[static_cast<std::size_t>(edge.second)].size()) - 1);
    }
    return counts;
}

TopologyCounts strand_topology_counts(const Settings& s) {
    if (s.strand_topology == "ring")
        return {s.n1, s.n1, s.n1};
    if (s.strand_topology == "star") {
        const StarArchitecture architecture = star_architecture(s);
        return graph_topology_counts(s.n1, architecture.edges);
    }
    if (s.strand_topology == "grafted") {
        const GraftedArchitecture architecture = grafted_architecture(s);
        return graph_topology_counts(s.n1, architecture.edges);
    }
    return {
        std::max(0, s.n1 - 1),
        std::max(0, s.n1 - 2),
        std::max(0, s.n1 - 3)
    };
}

void add_graph_topology(
    System& sys,
    int first_atom,
    int bead_count,
    const std::vector<std::pair<int, int>>& edges
) {
    const std::vector<std::vector<int>> adjacency =
        graph_adjacency(bead_count, edges);
    for (const auto& edge : edges)
        sys.bonds.push_back({
            static_cast<int>(sys.bonds.size()) + 1, 1,
            first_atom + edge.first, first_atom + edge.second
        });
    for (int center = 0; center < bead_count; ++center) {
        const auto& neighbors = adjacency[static_cast<std::size_t>(center)];
        for (std::size_t a = 0; a < neighbors.size(); ++a)
            for (std::size_t c = a + 1; c < neighbors.size(); ++c)
                sys.angles.push_back({
                    static_cast<int>(sys.angles.size()) + 1, 1,
                    first_atom + neighbors[a], first_atom + center,
                    first_atom + neighbors[c]
                });
    }
    for (const auto& edge : edges) {
        for (int a : adjacency[static_cast<std::size_t>(edge.first)]) {
            if (a == edge.second) continue;
            for (int d : adjacency[static_cast<std::size_t>(edge.second)]) {
                if (d == edge.first) continue;
                sys.dihedrals.push_back({
                    static_cast<int>(sys.dihedrals.size()) + 1, 1,
                    first_atom + a, first_atom + edge.first,
                    first_atom + edge.second, first_atom + d
                });
            }
        }
    }
}

std::vector<pdms_filler::Vec3> star_arm_directions(int arm_count) {
    using pdms_filler::Vec3;
    const double root_three = std::sqrt(3.0);
    const double branch_projection = arm_count == 8 ? 0.5 : -std::cos(
        pdms_filler::radians(pdms_filler::kBondAngleDegrees));
    const double radial = std::sqrt(
        1.0 - branch_projection*branch_projection);
    std::vector<Vec3> directions;
    if (arm_count == 3) {
        directions = {
            {1.0, 0.0, 0.0},
            {-0.5, 0.5*root_three, 0.0},
            {-0.5, -0.5*root_three, 0.0}
        };
    } else if (arm_count == 4) {
        directions = {
            {1.0, 1.0, 1.0}, {1.0, -1.0, -1.0},
            {-1.0, 1.0, -1.0}, {-1.0, -1.0, 1.0}
        };
    } else {
        for (int arm = 0; arm < 3; ++arm) {
            const double phi = 2.0*kPi*arm/3.0;
            directions.push_back({
                -branch_projection,
                radial*std::cos(phi), radial*std::sin(phi)
            });
        }
        if (arm_count == 8) {
            directions.push_back({0.0, 0.0, 1.0});
            directions.push_back({0.0, 0.0, -1.0});
        }
        const double phase = kPi/3.0;
        for (int arm = 0; arm < 3; ++arm) {
            const double phi = phase + 2.0*kPi*arm/3.0;
            directions.push_back({
                branch_projection,
                radial*std::cos(phi), radial*std::sin(phi)
            });
        }
    }
    for (Vec3& direction : directions)
        direction = pdms_filler::normalized(direction);
    return directions;
}

void add_linear_component(System& sys, int bead_count, int molecule_count, int component,
                          const Box& box, double bond_length, double spacing,
                          const std::set<int>& reactive_sites = {}) {
    if (bead_count == 0 && molecule_count != 0) throw std::runtime_error("A nonzero molecule count requires N > 0");
    if (molecule_count == 0) return;

    const double span = std::max(0, bead_count - 1) * bond_length;
    if (box.lx < span + 2.0*spacing)
        throw std::runtime_error("Lx is too small for a straight chain plus placement margins");
    if (box.ly < 2.0*spacing || box.lz < 2.0*spacing)
        throw std::runtime_error("Ly and Lz must each be at least twice the placement spacing");
    const int nx = std::max(1, static_cast<int>(box.lx / (span + spacing)));
    const int ny = std::max(1, static_cast<int>(box.ly / spacing));
    const int nz = std::max(1, static_cast<int>(box.lz / spacing));
    const long long capacity = 1LL * nx * ny * nz;
    if (capacity < molecule_count) {
        std::ostringstream msg;
        msg << "Placement grid capacity (" << capacity << ") is smaller than M" << component
            << " (" << molecule_count
            << "). Reduce molecule count or length, or lower the initial density.";
        throw std::runtime_error(msg.str());
    }

    for (int molecule_index = 0; molecule_index < molecule_count; ++molecule_index) {
        const int ix = molecule_index % nx;
        const int iy = (molecule_index / nx) % ny;
        const int iz = molecule_index / (nx * ny);
        const bool reverse = component != 1;
        double x = reverse ? box.lx / 2 - spacing - ix * (span + spacing)
                           : -box.lx / 2 + spacing + ix * (span + spacing);
        const double y = (reverse ? box.ly / 2 - spacing : -box.ly / 2 + spacing) + (reverse ? -1 : 1) * iy * spacing;
        double z = 0.0;
        if (component == 2) z = box.lz / 2 - spacing - iz * spacing;
        else z = -box.lz / 2 + spacing + iz * spacing;

        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 1; bead <= bead_count; ++bead) {
            int atom_type = 1;
            if (component == 1 && reactive_sites.count(bead)) atom_type = 2;
            if (component == 2 && reactive_sites.count(bead)) atom_type = 3;
            sys.atoms.push_back({static_cast<int>(sys.atoms.size()) + 1, molecule_id, atom_type, 0.0, x, y, z});
            x += reverse ? -bond_length : bond_length;
        }
        add_linear_topology(sys, first_atom, bead_count);
    }
}

void add_ring_strands(System& sys, const Settings& s, const Box& box,
                      const std::set<int>& reactive_sites) {
    const double radius =
        kBondLengthAngstrom / (2.0 * std::sin(kPi / s.n1));
    const double diameter = 2.0 * radius;
    const double spacing = kPlacementSpacing800KAngstrom;
    if (box.lx < diameter + 2.0*spacing ||
        box.ly < diameter + 2.0*spacing || box.lz < 2.0*spacing)
        throw std::runtime_error(
            "Simulation box is too small for the requested ring strands");

    const int nx = static_cast<int>(
        std::floor((box.lx - 2.0*spacing - diameter) /
                   (diameter + spacing))) + 1;
    const int ny = static_cast<int>(
        std::floor((box.ly - 2.0*spacing - diameter) /
                   (diameter + spacing))) + 1;
    const int nz = static_cast<int>(
        std::floor((box.lz - 2.0*spacing) / spacing)) + 1;
    const long long capacity = 1LL * nx * ny * nz;
    if (capacity < s.m1) {
        std::ostringstream message;
        message << "Ring placement capacity (" << capacity
                << ") is smaller than the strand count (" << s.m1
                << "). Reduce strand count or length, or lower initial density.";
        throw std::runtime_error(message.str());
    }

    for (int molecule_index = 0; molecule_index < s.m1; ++molecule_index) {
        const int ix = molecule_index % nx;
        const int iy = (molecule_index / nx) % ny;
        const int iz = molecule_index / (nx * ny);
        const double cx = -box.lx/2 + spacing + radius +
                          ix * (diameter + spacing);
        const double cy = -box.ly/2 + spacing + radius +
                          iy * (diameter + spacing);
        const double cz = -box.lz/2 + spacing + iz * spacing;
        const double phase = (molecule_index % 2) * kPi / s.n1;
        const int molecule_id = static_cast<int>(
            sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 1; bead <= s.n1; ++bead) {
            const double theta =
                2.0 * kPi * (bead - 1) / s.n1 + phase;
            const int atom_type = reactive_sites.count(bead) ? 2 : 1;
            sys.atoms.push_back({
                static_cast<int>(sys.atoms.size()) + 1, molecule_id,
                atom_type, 0.0, cx + radius*std::cos(theta),
                cy + radius*std::sin(theta), cz
            });
        }
        add_ring_topology(sys, first_atom, s.n1);
    }
}

bool internal_geometry_allowed(
    const std::vector<pdms_filler::Vec3>& positions,
    const std::vector<std::vector<int>>& adjacency
) {
    const int bead_count = static_cast<int>(positions.size());
    for (int source = 0; source < bead_count; ++source) {
        std::vector<int> distance(static_cast<std::size_t>(bead_count), -1);
        std::vector<int> queue;
        queue.reserve(static_cast<std::size_t>(bead_count));
        distance[static_cast<std::size_t>(source)] = 0;
        queue.push_back(source);
        for (std::size_t next = 0; next < queue.size(); ++next) {
            const int current = queue[next];
            for (int neighbor : adjacency[static_cast<std::size_t>(current)]) {
                if (distance[static_cast<std::size_t>(neighbor)] >= 0) continue;
                distance[static_cast<std::size_t>(neighbor)] =
                    distance[static_cast<std::size_t>(current)] + 1;
                queue.push_back(neighbor);
            }
        }
        for (int other = source + 1; other < bead_count; ++other) {
            const int path = distance[static_cast<std::size_t>(other)];
            if (path <= 2) continue;
            const double minimum = path == 3 ? 4.0 : 4.5;
            if (pdms_filler::norm2(
                    positions[static_cast<std::size_t>(source)] -
                    positions[static_cast<std::size_t>(other)]) <
                minimum * minimum)
                return false;
        }
    }
    return true;
}

std::vector<std::vector<int>> graph_distance_matrix(
    const std::vector<std::vector<int>>& adjacency
) {
    const int bead_count = static_cast<int>(adjacency.size());
    std::vector<std::vector<int>> distances(
        static_cast<std::size_t>(bead_count),
        std::vector<int>(static_cast<std::size_t>(bead_count), -1));
    for (int source = 0; source < bead_count; ++source) {
        std::vector<int> queue = {source};
        distances[static_cast<std::size_t>(source)]
                 [static_cast<std::size_t>(source)] = 0;
        for (std::size_t next = 0; next < queue.size(); ++next) {
            const int current = queue[next];
            for (int neighbor : adjacency[static_cast<std::size_t>(current)]) {
                int& distance = distances[static_cast<std::size_t>(source)]
                                         [static_cast<std::size_t>(neighbor)];
                if (distance >= 0) continue;
                distance = distances[static_cast<std::size_t>(source)]
                                    [static_cast<std::size_t>(current)] + 1;
                queue.push_back(neighbor);
            }
        }
    }
    return distances;
}

std::vector<pdms_filler::Vec3> build_local_star(
    const Settings& s,
    const StarArchitecture& architecture,
    const pdms_filler::Vec3& half_extent,
    std::mt19937& random
) {
    using pdms_filler::Vec3;
    const std::vector<Vec3> directions =
        star_arm_directions(s.star_arm_count);
    const std::vector<std::vector<int>> adjacency =
        graph_adjacency(s.n1, architecture.edges);
    const std::vector<std::vector<int>> distances =
        graph_distance_matrix(adjacency);
    std::uniform_real_distribution<double> twist(0.0, 2.0*kPi);
    int complete_candidates = 0;
    int internally_allowed_candidates = 0;
    int maximum_arms_built = 0;

    for (int attempt = 0; attempt < 500; ++attempt) {
        std::vector<Vec3> positions(static_cast<std::size_t>(s.n1));
        std::vector<bool> placed(static_cast<std::size_t>(s.n1), false);
        for (int center = 0; center < architecture.center_count; ++center) {
            positions[static_cast<std::size_t>(center)] = {
                (center - 0.5*(architecture.center_count - 1)) *
                    kBondLengthAngstrom,
                0.0,
                0.0
            };
            placed[static_cast<std::size_t>(center)] = true;
        }

        bool built = true;
        std::vector<int> arm_order;
        if (s.star_arm_count == 8)
            arm_order = {0, 1, 2, 5, 6, 7, 3, 4};
        else {
            for (int arm = 0; arm < s.star_arm_count; ++arm)
                arm_order.push_back(arm);
        }
        for (int order_index = 0;
             order_index < s.star_arm_count && built; ++order_index) {
            const int arm = arm_order[static_cast<std::size_t>(order_index)];
            const Vec3 direction = directions[static_cast<std::size_t>(arm)];
            const auto raw_basis = pdms_filler::perpendicular_basis(direction);
            const Vec3 center = positions[static_cast<std::size_t>(
                architecture.attachment_centers[static_cast<std::size_t>(arm)])];
            const int first = architecture.center_count +
                              arm * s.star_arm_length;
            bool arm_built = false;
            for (int arm_attempt = 0; arm_attempt < 500 && !arm_built;
                 ++arm_attempt) {
                const double phi = twist(random);
                const Vec3 first_basis =
                    std::cos(phi)*raw_basis.first +
                    std::sin(phi)*raw_basis.second;
                const Vec3 second_basis =
                    -std::sin(phi)*raw_basis.first +
                    std::cos(phi)*raw_basis.second;
                const int straight_prefix = s.star_arm_count == 8 ? 3 : 1;
                const double cone =
                    (s.star_arm_count == 6 ? 74.0 : 78.0)*kPi/180.0;
                const double turn =
                    (s.star_arm_count == 6 ? 60.0 : 45.0)*kPi/180.0;
                const double handedness = arm_attempt % 2 == 0 ? 1.0 : -1.0;
                Vec3 position = center;
                for (int bead = 1; bead <= s.star_arm_length; ++bead) {
                    Vec3 step_direction = direction;
                    if (bead > straight_prefix) {
                        const double azimuth =
                            handedness * (bead - straight_prefix - 1) * turn;
                        step_direction =
                            std::cos(cone)*direction +
                            std::sin(cone) *
                                (std::cos(azimuth)*first_basis +
                                 std::sin(azimuth)*second_basis);
                    }
                    position = position +
                               kBondLengthAngstrom*step_direction;
                    positions[static_cast<std::size_t>(first + bead - 1)] =
                        position;
                }

                bool allowed = true;
                for (int bead = first;
                     bead < first + s.star_arm_length && allowed; ++bead) {
                    for (int other = 0; other < s.n1; ++other) {
                        if (other == bead) continue;
                        const bool current_arm =
                            other >= first &&
                            other < first + s.star_arm_length;
                        if (!placed[static_cast<std::size_t>(other)] &&
                            !current_arm)
                            continue;
                        const int path =
                            distances[static_cast<std::size_t>(bead)]
                                     [static_cast<std::size_t>(other)];
                        if (path <= 2) continue;
                        const double minimum = path == 3 ? 4.0 : 4.5;
                        if (pdms_filler::norm2(
                                positions[static_cast<std::size_t>(bead)] -
                                positions[static_cast<std::size_t>(other)]) <
                            minimum*minimum) {
                            allowed = false;
                            break;
                        }
                    }
                }
                arm_built = allowed;
            }
            if (!arm_built) {
                built = false;
                break;
            }
            for (int bead = first; bead < first + s.star_arm_length; ++bead)
                placed[static_cast<std::size_t>(bead)] = true;
            maximum_arms_built = std::max(
                maximum_arms_built, order_index + 1);
        }
        if (!built) continue;
        ++complete_candidates;
        if (!internal_geometry_allowed(positions, adjacency))
            continue;
        ++internally_allowed_candidates;

        Vec3 lower = positions.front();
        Vec3 upper = positions.front();
        for (const Vec3& position : positions) {
            lower.x = std::min(lower.x, position.x);
            lower.y = std::min(lower.y, position.y);
            lower.z = std::min(lower.z, position.z);
            upper.x = std::max(upper.x, position.x);
            upper.y = std::max(upper.y, position.y);
            upper.z = std::max(upper.z, position.z);
        }
        const Vec3 midpoint = 0.5*(lower + upper);
        for (Vec3& position : positions) position = position - midpoint;
        if (0.5*(upper.x - lower.x) > half_extent.x ||
            0.5*(upper.y - lower.y) > half_extent.y ||
            0.5*(upper.z - lower.z) > half_extent.z)
            continue;
        return positions;
    }
    std::ostringstream message;
    message << "Could not build a compact, self-avoiding star conformation "
            << "for the available placement cell (complete candidates: "
            << complete_candidates << ", internally allowed: "
            << internally_allowed_candidates << ", half extents: "
            << half_extent.x << ", " << half_extent.y << ", "
            << half_extent.z << " A, maximum arms built: "
            << maximum_arms_built
            << "). Reduce star count or arm length, or "
               "lower the initial density.";
    throw std::runtime_error(message.str());
}

std::vector<pdms_filler::Vec3> compact_grafted_backbone(
    int bead_count,
    int beads_per_turn
) {
    using pdms_filler::Vec3;
    constexpr double pitch = 4.8;
    const double z_step = pitch / beads_per_turn;
    const double theta_step = 2.0*kPi / beads_per_turn;
    const double lateral_step = std::sqrt(
        kBondLengthAngstrom*kBondLengthAngstrom - z_step*z_step);
    const double radius = lateral_step /
        (2.0*std::sin(0.5*theta_step));
    std::vector<Vec3> backbone(static_cast<std::size_t>(bead_count));
    for (int bead = 0; bead < bead_count; ++bead) {
        const double theta = bead*theta_step;
        backbone[static_cast<std::size_t>(bead)] = {
            radius*std::cos(theta), radius*std::sin(theta), bead*z_step
        };
    }
    const double z_midpoint = 0.5*(bead_count - 1)*z_step;
    for (Vec3& position : backbone) position.z -= z_midpoint;
    return backbone;
}

bool graft_candidate_allowed(
    const std::vector<pdms_filler::Vec3>& positions,
    const std::vector<bool>& placed,
    int first,
    int count,
    const std::vector<std::vector<int>>& distances
) {
    for (int bead = first; bead < first + count; ++bead) {
        for (int other = 0; other < static_cast<int>(positions.size());
             ++other) {
            if (other == bead) continue;
            const bool current_side =
                other >= first && other < first + count;
            if (!placed[static_cast<std::size_t>(other)] && !current_side)
                continue;
            const int path = distances[static_cast<std::size_t>(bead)]
                                      [static_cast<std::size_t>(other)];
            if (path <= 2) continue;
            const double minimum = path == 3 ? 4.0 : 4.5;
            if (pdms_filler::norm2(
                    positions[static_cast<std::size_t>(bead)] -
                    positions[static_cast<std::size_t>(other)]) <
                minimum*minimum)
                return false;
        }
    }
    return true;
}

std::vector<pdms_filler::Vec3> build_local_grafted(
    const Settings& s,
    const GraftedArchitecture& architecture,
    const pdms_filler::Vec3& half_extent,
    std::mt19937& random
) {
    using pdms_filler::Vec3;
    const std::vector<std::vector<int>> adjacency =
        graph_adjacency(s.n1, architecture.edges);
    const std::vector<std::vector<int>> distances =
        graph_distance_matrix(adjacency);
    std::uniform_real_distribution<double> phase(0.0, 2.0*kPi);
    int complete_candidates = 0;
    int internally_allowed_candidates = 0;
    int maximum_sides_built = 0;

    for (int attempt = 0; attempt < 500; ++attempt) {
        std::vector<Vec3> positions(static_cast<std::size_t>(s.n1));
        std::vector<bool> placed(static_cast<std::size_t>(s.n1), false);
        const std::vector<Vec3> backbone =
            compact_grafted_backbone(
                s.graft_backbone_length, s.graft_spacing == 0 ? 14 : 16);
        for (int bead = 0; bead < s.graft_backbone_length; ++bead) {
            positions[static_cast<std::size_t>(bead)] =
                backbone[static_cast<std::size_t>(bead)];
            placed[static_cast<std::size_t>(bead)] = true;
        }

        bool built = true;
        for (int side = 0; side < s.graft_side_chain_count && built; ++side) {
            const int graft =
                architecture.graft_sites[static_cast<std::size_t>(side)];
            const Vec3 center = backbone[static_cast<std::size_t>(graft)];
            const Vec3 direction = pdms_filler::normalized(
                Vec3{center.x, center.y, 0.0});
            const auto raw_basis = pdms_filler::perpendicular_basis(direction);
            const int first = s.graft_backbone_length +
                              side*s.graft_side_chain_length;
            bool side_built = false;
            for (int side_attempt = 0;
                 side_attempt < 500 && !side_built; ++side_attempt) {
                const double phi = phase(random);
                const Vec3 first_basis =
                    std::cos(phi)*raw_basis.first +
                    std::sin(phi)*raw_basis.second;
                const Vec3 second_basis =
                    -std::sin(phi)*raw_basis.first +
                    std::cos(phi)*raw_basis.second;
                const bool dense = s.graft_spacing == 0;
                const int straight_prefix = dense
                    ? s.graft_side_chain_length : 1;
                const double cone = 78.0*kPi/180.0;
                const double turn = 45.0*kPi/180.0;
                const double handedness = side_attempt % 2 == 0 ? 1.0 : -1.0;
                Vec3 position = center;
                for (int bead = 1; bead <= s.graft_side_chain_length; ++bead) {
                    Vec3 step_direction = direction;
                    if (bead > straight_prefix) {
                        const double azimuth = handedness *
                            (bead - straight_prefix - 1)*turn;
                        step_direction =
                            std::cos(cone)*direction +
                            std::sin(cone)*
                                (std::cos(azimuth)*first_basis +
                                 std::sin(azimuth)*second_basis);
                    }
                    position = position +
                               kBondLengthAngstrom*step_direction;
                    positions[static_cast<std::size_t>(first + bead - 1)] =
                        position;
                }
                side_built = graft_candidate_allowed(
                    positions, placed, first, s.graft_side_chain_length,
                    distances);
            }
            if (!side_built) {
                built = false;
                break;
            }
            for (int bead = first;
                 bead < first + s.graft_side_chain_length; ++bead)
                placed[static_cast<std::size_t>(bead)] = true;
            maximum_sides_built = std::max(maximum_sides_built, side + 1);
        }
        if (!built) continue;
        ++complete_candidates;
        if (!internal_geometry_allowed(positions, adjacency)) continue;
        ++internally_allowed_candidates;

        Vec3 lower = positions.front();
        Vec3 upper = positions.front();
        for (const Vec3& position : positions) {
            lower.x = std::min(lower.x, position.x);
            lower.y = std::min(lower.y, position.y);
            lower.z = std::min(lower.z, position.z);
            upper.x = std::max(upper.x, position.x);
            upper.y = std::max(upper.y, position.y);
            upper.z = std::max(upper.z, position.z);
        }
        const Vec3 midpoint = 0.5*(lower + upper);
        for (Vec3& position : positions) position = position - midpoint;
        if (0.5*(upper.x - lower.x) > half_extent.x ||
            0.5*(upper.y - lower.y) > half_extent.y ||
            0.5*(upper.z - lower.z) > half_extent.z)
            continue;
        return positions;
    }

    std::ostringstream message;
    message << "Could not build a compact, self-avoiding grafted conformation "
            << "for the available placement cell (complete candidates: "
            << complete_candidates << ", internally allowed: "
            << internally_allowed_candidates << ", half extents: "
            << half_extent.x << ", " << half_extent.y << ", "
            << half_extent.z << " A, maximum side chains built: "
            << maximum_sides_built
            << "). Reduce molecule count or chain lengths, increase graft "
               "spacing, or lower the initial density.";
    throw std::runtime_error(message.str());
}

struct BranchedPlacementGrid {
    int nx = 1;
    int ny = 1;
    int nz = 1;
    double cell_x = 0.0;
    double cell_y = 0.0;
    double cell_z = 0.0;
    double z_lower = 0.0;
};

BranchedPlacementGrid branched_placement_grid(
    const Settings& s,
    const Box& box,
    double z_lower,
    double z_upper
) {
    const double z_length = z_upper - z_lower;
    if (z_length <= kPlacementSpacing800KAngstrom)
        throw std::runtime_error(
            "Insufficient z space to place branched strands below crosslinkers");
    const int base = std::max(
        1, static_cast<int>(std::ceil(std::cbrt(s.m1))));
    const int maximum = 4 * base;
    double best_score = -1.0;
    long long best_capacity = 0;
    BranchedPlacementGrid best;
    for (int nx = 1; nx <= maximum; ++nx) {
        for (int ny = 1; ny <= maximum; ++ny) {
            const long long plane = 1LL*nx*ny;
            const int nz = static_cast<int>((s.m1 + plane - 1) / plane);
            if (nz < 1 || nz > maximum) continue;
            const double cell_x = box.lx / nx;
            const double cell_y = box.ly / ny;
            const double cell_z = z_length / nz;
            const double score = std::min({cell_x, cell_y, cell_z});
            const long long capacity = plane*nz;
            if (score > best_score + 1.0e-12 ||
                (std::fabs(score - best_score) <= 1.0e-12 &&
                 (best_capacity == 0 || capacity < best_capacity))) {
                best_score = score;
                best_capacity = capacity;
                best = {nx, ny, nz, cell_x, cell_y, cell_z, z_lower};
            }
        }
    }
    if (best_capacity < s.m1)
        throw std::runtime_error("Could not construct a branched-strand placement grid");
    return best;
}

double branched_placement_upper_z(const Settings& s, const Box& box) {
    if (s.m2 == 0) return box.lz/2;
    const double span = std::max(0, s.n2 - 1) * kBondLengthAngstrom;
    const int nx = std::max(
        1, static_cast<int>(box.lx /
                            (span + kPlacementSpacing800KAngstrom)));
    const int ny = std::max(
        1, static_cast<int>(box.ly / kPlacementSpacing800KAngstrom));
    const long long per_layer = 1LL*nx*ny;
    const long long layers = (s.m2 + per_layer - 1) / per_layer;
    const double lowest_crosslinker_z =
        box.lz/2 - kPlacementSpacing800KAngstrom -
        (layers - 1)*kPlacementSpacing800KAngstrom;
    return lowest_crosslinker_z - kPlacementSpacing800KAngstrom;
}

void add_star_strands(
    System& sys,
    const Settings& s,
    const Box& box,
    const std::set<int>& reactive_sites
) {
    const StarArchitecture architecture = star_architecture(s);
    const double z_lower = -box.lz/2;
    const double z_upper = branched_placement_upper_z(s, box);
    const BranchedPlacementGrid grid =
        branched_placement_grid(s, box, z_lower, z_upper);
    const double spacing = kPlacementSpacing800KAngstrom;
    const pdms_filler::Vec3 half_extent{
        0.5*(grid.cell_x - spacing),
        0.5*(grid.cell_y - spacing),
        0.5*(grid.cell_z - spacing)
    };
    if (half_extent.x <= 0.0 || half_extent.y <= 0.0 ||
        half_extent.z <= 0.0)
        throw std::runtime_error(
            "Star placement cells are smaller than the 7.5 A model spacing");
    std::mt19937 random(s.seed);
    const std::vector<pdms_filler::Vec3> local =
        build_local_star(s, architecture, half_extent, random);

    for (int molecule_index = 0; molecule_index < s.m1; ++molecule_index) {
        const int ix = molecule_index % grid.nx;
        const int iy = (molecule_index / grid.nx) % grid.ny;
        const int iz = molecule_index / (grid.nx * grid.ny);
        const pdms_filler::Vec3 translation{
            -box.lx/2 + (ix + 0.5)*grid.cell_x,
            -box.ly/2 + (iy + 0.5)*grid.cell_y,
            grid.z_lower + (iz + 0.5)*grid.cell_z
        };
        const int molecule_id = static_cast<int>(
            sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 0; bead < s.n1; ++bead) {
            const pdms_filler::Vec3 position =
                local[static_cast<std::size_t>(bead)] + translation;
            sys.atoms.push_back({
                static_cast<int>(sys.atoms.size()) + 1, molecule_id,
                reactive_sites.count(bead + 1) ? 2 : 1, 0.0,
                position.x, position.y, position.z
            });
        }
        add_graph_topology(
            sys, first_atom, s.n1, architecture.edges);
    }
}

void add_grafted_strands(
    System& sys,
    const Settings& s,
    const Box& box,
    const std::set<int>& reactive_sites
) {
    const GraftedArchitecture architecture = grafted_architecture(s);
    const double z_lower = -box.lz/2;
    const double z_upper = branched_placement_upper_z(s, box);
    const BranchedPlacementGrid grid =
        branched_placement_grid(s, box, z_lower, z_upper);
    const double spacing = kPlacementSpacing800KAngstrom;
    const pdms_filler::Vec3 half_extent{
        0.5*(grid.cell_x - spacing),
        0.5*(grid.cell_y - spacing),
        0.5*(grid.cell_z - spacing)
    };
    if (half_extent.x <= 0.0 || half_extent.y <= 0.0 ||
        half_extent.z <= 0.0)
        throw std::runtime_error(
            "Grafted placement cells are smaller than the 7.5 A model spacing");
    std::mt19937 random(s.seed);
    const std::vector<pdms_filler::Vec3> local =
        build_local_grafted(s, architecture, half_extent, random);

    for (int molecule_index = 0; molecule_index < s.m1; ++molecule_index) {
        const int ix = molecule_index % grid.nx;
        const int iy = (molecule_index / grid.nx) % grid.ny;
        const int iz = molecule_index / (grid.nx * grid.ny);
        const pdms_filler::Vec3 translation{
            -box.lx/2 + (ix + 0.5)*grid.cell_x,
            -box.ly/2 + (iy + 0.5)*grid.cell_y,
            grid.z_lower + (iz + 0.5)*grid.cell_z
        };
        const int molecule_id = static_cast<int>(
            sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 0; bead < s.n1; ++bead) {
            const pdms_filler::Vec3 position =
                local[static_cast<std::size_t>(bead)] + translation;
            sys.atoms.push_back({
                static_cast<int>(sys.atoms.size()) + 1, molecule_id,
                reactive_sites.count(bead + 1) ? 2 : 1, 0.0,
                position.x, position.y, position.z
            });
        }
        add_graph_topology(sys, first_atom, s.n1, architecture.edges);
    }
}

double minimum_image(double delta, double length) {
    return delta - std::round(delta / length) * length;
}

bool overlaps_existing(const std::vector<pdms_filler::Vec3>& candidate,
                       const System& sys, const Box& box,
                       double minimum_separation, bool periodic_z) {
    const double minimum_squared = minimum_separation * minimum_separation;
    for (const pdms_filler::Vec3& position : candidate) {
        for (const Atom& atom : sys.atoms) {
            double dx = minimum_image(position.x - atom.x, box.lx);
            double dy = minimum_image(position.y - atom.y, box.ly);
            double dz = position.z - atom.z;
            if (periodic_z) dz = minimum_image(dz, box.lz);
            if (dx*dx + dy*dy + dz*dz < minimum_squared)
                return true;
        }
    }
    return false;
}

void add_star_moderators(System& sys, const Settings& s, const Box& box,
                         std::mt19937& rng) {
    const double margin = kBondLengthAngstrom;
    if (box.lx <= 2.0*margin || box.ly <= 2.0*margin)
        throw std::runtime_error("Lateral box dimensions are too small for a star moderator");
    std::uniform_real_distribution<double> xdist(-box.lx/2 + margin, box.lx/2 - margin);
    std::uniform_real_distribution<double> ydist(-box.ly/2 + margin, box.ly/2 - margin);
    std::uniform_real_distribution<double> zdist(
        kModeratorZLowerFraction * box.lz,
        kModeratorZUpperFraction * box.lz);
    for (int i = 0; i < s.m3; ++i) {
        std::vector<pdms_filler::Vec3> candidate;
        bool accepted = false;
        for (int attempt = 0; attempt < 10000 && !accepted; ++attempt) {
            const double x = xdist(rng);
            const double y = ydist(rng);
            const double z = zdist(rng);
            const double b = kBondLengthAngstrom;
            candidate = {
                {x,     y,     z},
                {x + b, y,     z},
                {x - b, y,     z},
                {x,     y + b, z},
                {x,     y - b, z}
            };
            accepted = !overlaps_existing(
                candidate, sys, box, s.filler_minimum_separation,
                s.thickness <= 0.0);
        }
        if (!accepted) {
            std::ostringstream message;
            message << "Could not place moderator " << i + 1
                    << " without overlap in the upper-middle z region. "
                       "Reduce the minimum separation, component loading, "
                       "or initial density.";
            throw std::runtime_error(message.str());
        }

        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int center = static_cast<int>(sys.atoms.size()) + 1;
        // Preserve the legacy force-field mapping: ordinary center, type-2 arms.
        for (std::size_t bead = 0; bead < candidate.size(); ++bead) {
            const pdms_filler::Vec3& position = candidate[bead];
            sys.atoms.push_back({
                center + static_cast<int>(bead),
                molecule_id,
                bead == 0 ? 1 : 2,
                0.0,
                position.x,
                position.y,
                position.z
            });
        }
        for (int arm = 1; arm <= 4; ++arm)
            sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 1, center, center + arm});
        for (int a = 1; a <= 4; ++a)
            for (int c = a + 1; c <= 4; ++c)
                sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1, center + a, center, center + c});
    }
}

void add_pdms_filler(System& sys, const Settings& s, const Box& box) {
    if (s.m4 == 0) return;
    std::vector<pdms_filler::Vec3> existing;
    existing.reserve(sys.atoms.size());
    for (const Atom& atom : sys.atoms)
        existing.push_back({atom.x, atom.y, atom.z});

    pdms_filler::Settings filler;
    filler.length = s.n4;
    filler.chains = s.m4;
    filler.seed = s.filler_seed;
    filler.minimum_separation = s.filler_minimum_separation;
    filler.z_lower_fraction = kFillerZLowerFraction;
    filler.z_upper_fraction = kFillerZUpperFraction;
    const pdms_filler::Box filler_box{
        box.lx, box.ly, box.lz, s.thickness <= 0.0
    };
    const pdms_filler::Component generated =
        pdms_filler::generate(filler, filler_box, existing);

    const int atom_offset = static_cast<int>(sys.atoms.size());
    const int molecule_offset =
        sys.atoms.empty() ? 0 : sys.atoms.back().molecule;
    for (const pdms_filler::Atom& atom : generated.atoms) {
        sys.atoms.push_back({
            static_cast<int>(sys.atoms.size()) + 1,
            molecule_offset + atom.molecule,
            atom.type,
            0.0,
            atom.position.x,
            atom.position.y,
            atom.position.z
        });
    }
    for (const pdms_filler::Bond& bond : generated.bonds)
        sys.bonds.push_back({
            static_cast<int>(sys.bonds.size()) + 1, bond.type,
            atom_offset + bond.a, atom_offset + bond.b
        });
    for (const pdms_filler::Angle& angle : generated.angles)
        sys.angles.push_back({
            static_cast<int>(sys.angles.size()) + 1, angle.type,
            atom_offset + angle.a, atom_offset + angle.b,
            atom_offset + angle.c
        });
    for (const pdms_filler::Dihedral& dihedral : generated.dihedrals)
        sys.dihedrals.push_back({
            static_cast<int>(sys.dihedrals.size()) + 1, dihedral.type,
            atom_offset + dihedral.a, atom_offset + dihedral.b,
            atom_offset + dihedral.c, atom_offset + dihedral.d
        });
}

System build_system(const Settings& s, const Box& box) {
    System sys;
    const long long expected_atoms =
        1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n3*s.m3 + filler_beads(s);
    sys.atoms.reserve(static_cast<size_t>(expected_atoms));
    const std::set<int> strand_reactive_sites = strand_sites(s);
    std::cerr << "Strand type-2 sites (" << s.strand_topology << ", "
              << strand_reactive_distribution_name(s) << "):";
    for (int site : strand_reactive_sites) std::cerr << ' ' << site;
    std::cerr << '\n';
    if (s.strand_topology == "ring")
        add_ring_strands(sys, s, box, strand_reactive_sites);
    else if (s.strand_topology == "star")
        add_star_strands(sys, s, box, strand_reactive_sites);
    else if (s.strand_topology == "grafted")
        add_grafted_strands(sys, s, box, strand_reactive_sites);
    else
        add_linear_component(
            sys, s.n1, s.m1, 1, box, kBondLengthAngstrom,
            kPlacementSpacing800KAngstrom, strand_reactive_sites);
    const std::set<int> reactive_sites = crosslinker_sites(s);
    std::cerr << "Cross-linker type-3 sites (" << s.crosslink_distribution << "):";
    for (int site : reactive_sites) std::cerr << ' ' << site;
    std::cerr << '\n';
    add_linear_component(
        sys, s.n2, s.m2, 2, box, kBondLengthAngstrom,
        kPlacementSpacing800KAngstrom, reactive_sites);
    std::mt19937 rng(s.seed);
    add_star_moderators(sys, s, box, rng);
    add_pdms_filler(sys, s, box);
    if (static_cast<long long>(sys.atoms.size()) != expected_atoms)
        throw std::logic_error("Internal error: generated atom count does not match requested composition");
    return sys;
}

void write_data(const Settings& s, const System& sys, const Box& box,
                const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Cannot open output file: " + output_path);
    out << "LAMMPS data file for a coarse-grained PDMS elastomer\n\n"
        << sys.atoms.size() << " atoms\n3 atom types\n"
        << sys.bonds.size() << " bonds\n2 bond types\n"
        << sys.angles.size() << " angles\n1 angle types\n"
        << sys.dihedrals.size() << " dihedrals\n1 dihedral types\n\n"
        << std::fixed << std::setprecision(6)
        << -box.lx/2 << ' ' << box.lx/2 << " xlo xhi\n"
        << -box.ly/2 << ' ' << box.ly/2 << " ylo yhi\n"
        << -box.lz/2 << ' ' << box.lz/2 << " zlo zhi\n\nMasses\n\n";
    out << "1 " << s.bead_mass << '\n'
        << "2 " << s.bead_mass << '\n'
        << "3 " << s.bead_mass << '\n';
    out << "\nAtoms # full\n\n" << std::setprecision(8);
    for (const auto& a : sys.atoms)
        out << a.id << ' ' << a.molecule << ' ' << a.type << ' ' << a.charge << ' '
            << a.x << ' ' << a.y << ' ' << a.z << '\n';
    if (!sys.bonds.empty()) {
        out << "\nBonds\n\n";
        for (const auto& b : sys.bonds) out << b.id << ' ' << b.type << ' ' << b.a << ' ' << b.b << '\n';
    }
    if (!sys.angles.empty()) {
        out << "\nAngles\n\n";
        for (const auto& a : sys.angles) out << a.id << ' ' << a.type << ' ' << a.a << ' ' << a.b << ' ' << a.c << '\n';
    }
    if (!sys.dihedrals.empty()) {
        out << "\nDihedrals\n\n";
        for (const auto& d : sys.dihedrals) out << d.id << ' ' << d.type << ' ' << d.a << ' ' << d.b << ' ' << d.c << ' ' << d.d << '\n';
    }
    if (!out) throw std::runtime_error("Failed while writing output file: " + output_path);
}

struct OutputFiles {
    std::string directory;
    std::string data;
    std::string data_basename;
    std::string input;
    std::string input_basename;
    std::string submit;
    std::string submit_basename;
    std::string info;
    std::string info_basename;
    std::string case_name;
};

struct LjParameters {
    double epsilon;
    double sigma;
    double cutoff;
};

OutputFiles output_files(const Settings& s) {
    OutputFiles files;
    const std::filesystem::path requested_data(s.output);
    files.data_basename = requested_data.filename().string();
    if (files.data_basename.empty())
        throw std::runtime_error("Output filename cannot end with a directory separator");
    files.case_name = files.data_basename.rfind("data.", 0) == 0
        ? files.data_basename.substr(5) : files.data_basename;
    if (files.case_name.empty() || files.case_name == "." || files.case_name == "..")
        throw std::runtime_error("Cannot derive a case name from the output filename");

    const std::filesystem::path directory =
        requested_data.parent_path() / files.case_name;
    files.directory = directory.string();
    files.input_basename = "in." + files.case_name;
    files.submit_basename = "submit." + files.case_name + ".sh";
    files.info_basename = files.case_name + ".info";
    files.data = (directory / files.data_basename).string();
    files.input = (directory / files.input_basename).string();
    files.submit = (directory / files.submit_basename).string();
    files.info = (directory / files.info_basename).string();
    return files;
}

void create_output_directory(const OutputFiles& files) {
    std::error_code error;
    std::filesystem::create_directories(files.directory, error);
    if (error)
        throw std::runtime_error(
            "Cannot create output directory " + files.directory + ": " +
            error.message());
}

LjParameters lj_parameters(double temperature) {
    const double a = temperature - 186.04682;
    const double b = 0.00758 * a;
    const double c = 1.0 + std::exp(b);
    const double epsilon = (4.77795 / c + 1.47169) * 0.350646;
    const double sigma = ((7.86548e-05) * temperature + 1.27856) * 4.95013;
    return {epsilon, sigma, sigma * std::pow(2.0, 1.0 / 6.0)};
}

std::string sanitize_job_name(const std::string& value) {
    std::string result;
    for (char c : value) {
        const unsigned char uc = static_cast<unsigned char>(c);
        result.push_back(std::isalnum(uc) || c == '_' || c == '-' ? c : '_');
    }
    if (result.size() > 100) result.resize(100);
    return result.empty() ? "PDMS" : result;
}

std::string shell_single_quote(const std::string& value) {
    std::string result = "'";
    for (char c : value) {
        if (c == '\'') result += "'\\''";
        else result.push_back(c);
    }
    result += "'";
    return result;
}

std::string json_escape(const std::string& value) {
    std::ostringstream escaped;
    for (unsigned char c : value) {
        switch (c) {
            case '"': escaped << "\\\""; break;
            case '\\': escaped << "\\\\"; break;
            case '\b': escaped << "\\b"; break;
            case '\f': escaped << "\\f"; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default:
                if (c < 0x20)
                    escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                            << static_cast<int>(c) << std::dec << std::setfill(' ');
                else
                    escaped << static_cast<char>(c);
        }
    }
    return escaped.str();
}

void write_lammps_input(const Settings& s, const OutputFiles& files) {
    std::ofstream out(files.input);
    if (!out) throw std::runtime_error("Cannot open LAMMPS input file: " + files.input);

    const bool film = s.thickness > 0.0;
    const double compression_scale = film
        ? std::sqrt(s.density / s.target_density)
        : std::cbrt(s.density / s.target_density);
    const double initial_to_curing_scale = film
        ? std::sqrt(s.density / kControlledCuringDensity)
        : std::cbrt(s.density / kControlledCuringDensity);
    const double curing_to_target_scale = film
        ? std::sqrt(kControlledCuringDensity / s.target_density)
        : std::cbrt(kControlledCuringDensity / s.target_density);
    const LjParameters hot = lj_parameters(800.0);
    const LjParameters cold = lj_parameters(300.0);
    const double hot_global_cutoff =
        pdms_filler::maximum_repulsive_cutoff(800.0);
    const std::string suffix = files.case_name;
    const bool controlled_conversion = conversion_control_enabled(s);
    const long long maximum_new_bonds = stoichiometric_maximum_bonds(s);
    const long long requested_new_bonds = target_new_bonds(s);

    out << "# Generated by the generic PDMS elastomer generator\n"
        << "# Geometry: " << (film ? "film with fixed Lz" : "periodic bulk") << "\n\n"
        << "units           real\n"
        << "boundary        p p " << (film ? "f" : "p") << "\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     harmonic\n"
        << "dihedral_style  nharmonic\n"
        << "special_bonds   lj 0 0 0.5\n"
        << std::fixed << std::setprecision(8)
        << "pair_style      lj/cut " << hot_global_cutoff << "\n"
        << "comm_modify     cutoff 25\n"
        << "read_data       " << files.data_basename
        << " extra/bond/per/atom 4 extra/angle/per/atom 10"
        << " extra/dihedral/per/atom 10 extra/special/per/atom 30\n\n"
        << "mass            1 " << s.bead_mass << "\n"
        << "mass            2 " << s.bead_mass << "\n"
        << "mass            3 " << s.bead_mass << "\n\n"
        << "bond_coeff      1 115.4086 " << kBondLengthAngstrom << "\n"
        << "bond_coeff      2 115.4086 7.235\n"
        << "angle_coeff     1 64.62431 111.623\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n\n"
        << "# Explicit 800 K repulsive PDMS pair matrix\n";
    pdms_filler::write_pair_matrix(out, 800.0, true);
    out << "\n"
        << "neighbor        2 bin\n"
        << "neigh_modify    delay 0 every 1 check yes\n"
        << "timestep        5\n"
        << "thermo          1000\n"
        << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n\n";

    if (film) {
        out << "# Separate repulsive fixes for the lower and upper z box edges.\n"
            << "fix             zlo_wall all wall/lj126 zlo EDGE "
            << hot.epsilon << ' ' << hot.sigma << ' ' << hot.cutoff << " units box\n"
            << "fix             zhi_wall all wall/lj126 zhi EDGE "
            << hot.epsilon << ' ' << hot.sigma << ' ' << hot.cutoff << " units box\n"
            << "fix_modify      zlo_wall energy yes\n"
            << "fix_modify      zhi_wall energy yes\n\n";
    }

    out << "# Relax initial overlaps before assigning 800 K velocities.\n"
        << "min_style       cg\n"
        << "min_modify      dmax 0.1 line backtrack\n"
        << "minimize        1.0e-4 1.0e-4 1000 10000\n";
    if (film) {
        out << "fix_modify      zlo_wall energy no\n"
            << "fix_modify      zhi_wall energy no\n";
    }
    out << "reset_timestep  0\n"
        << "neigh_modify    delay 5 every 1 check yes\n\n"
        << "restart         100000 restart." << suffix << ".1 restart." << suffix << ".2\n"
        << "velocity        all create 800.0 " << s.seed
        << " mom yes rot yes dist gaussian\n\n"
        << "# 800 K relaxation\n"
        << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
        << "run             1000000\n"
        << "write_data      data." << suffix << ".rep_800 nocoeff\n\n"
        << "unfix           integrate\n\n";

    if (controlled_conversion) {
        out << "# Compress to the 0.5 g/cm^3 curing density without reactions\n"
            << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
            << "fix             compress_to_cure all deform 1 x scale "
            << initial_to_curing_scale << " y scale " << initial_to_curing_scale;
        if (!film) out << " z scale " << initial_to_curing_scale;
        out << " units box\n"
            << "run             " << kControlledCompressionStageSteps << "\n"
            << "unfix           compress_to_cure\n\n"
            << "# Cure at 0.5 g/cm^3 until conversion is reached or 5M steps elapse\n"
            << "fix             xlink all bond/create 1 2 3 " << hot.cutoff
            << " 2 iparam 1 1 jparam 1 1 prob "
            << kControlledCrosslinkProbability << " 348154\n"
            << "variable        target_new_bonds equal " << requested_new_bonds << "\n"
            << "variable        stoichiometric_maximum_bonds equal "
            << maximum_new_bonds << "\n"
            << "variable        created_new_bonds equal f_xlink[2]\n"
            << "variable        conversion_percent equal 100.0*v_created_new_bonds/v_stoichiometric_maximum_bonds\n"
            << "thermo_style    custom step temp density etotal epair ebond eangle"
            << " edihed f_xlink[1] f_xlink[2] v_conversion_percent bonds\n"
            << "fix             conversion_halt all halt 1 v_created_new_bonds >= "
            << requested_new_bonds << " error continue\n"
            << "run             " << kConversionControlledHoldMaximumSteps << "\n"
            << "unfix           conversion_halt\n"
            << "print           \"Conversion-controlled curing ended with "
               "$(f_xlink[2]:%.0f) new bonds; target "
            << requested_new_bonds << " of " << maximum_new_bonds << "\"\n"
            << "unfix           xlink\n"
            << "variable        conversion_percent delete\n"
            << "variable        created_new_bonds delete\n"
            << "variable        stoichiometric_maximum_bonds delete\n"
            << "variable        target_new_bonds delete\n\n"
            << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
            << " etotal epair ebond eangle edihed\n\n"
            << "# Compress the cured network from 0.5 g/cm^3 to the target density\n"
            << "fix             compress_to_target all deform 1 x scale "
            << curing_to_target_scale << " y scale " << curing_to_target_scale;
        if (!film) out << " z scale " << curing_to_target_scale;
        out << " units box\n"
            << "run             " << kControlledCompressionStageSteps << "\n"
            << "unfix           compress_to_target\n\n"
            << "# Equilibrate the cured network at 800 K without further reactions\n"
            << "run             " << kControlledPostCureEquilibrationSteps << "\n"
            << "unfix           integrate\n"
            << "reset_timestep  0\n";
    } else if (film) {
        out << "# Lateral compression at fixed box Lz and nominal wall-free film thickness\n"
            << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
            << "fix             xlink all bond/create 1 2 3 " << hot.cutoff
            << " 2 iparam 1 1 jparam 1 1 prob "
            << s.crosslink_probability << " 348154\n"
            << "thermo_style    custom step temp density etotal epair ebond eangle"
            << " edihed f_xlink[1] f_xlink[2] bonds\n"
            << "fix             compress all deform 1 x scale " << compression_scale
            << " y scale " << compression_scale << " units box\n"
            << "run             1000000\n"
            << "unfix           compress\n\n"
            << "# Relax at the compressed dimensions with continued crosslinking\n"
            << "run             1000000\n\n"
            << "# Continue crosslinking at fixed dimensions\n"
            << "run             2000000\n"
            << "unfix           integrate\n"
            << "unfix           xlink\n";
    } else {
        out << "# Isotropic compression with crosslinking at 800 K\n"
            << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
            << "fix             xlink all bond/create 1 2 3 " << hot.cutoff
            << " 2 iparam 1 1 jparam 1 1 prob "
            << s.crosslink_probability << " 348154\n"
            << "thermo_style    custom step temp density etotal epair ebond eangle"
            << " edihed f_xlink[1] f_xlink[2] bonds\n"
            << "fix             compress all deform 1 x scale " << compression_scale
            << " y scale " << compression_scale << " z scale " << compression_scale
            << " units box\n"
            << "run             1000000\n"
            << "unfix           compress\n\n"
            << "# Relax at the compressed dimensions with continued crosslinking\n"
            << "run             1000000\n\n"
            << "# Continue crosslinking at fixed dimensions\n"
            << "run             2000000\n"
            << "unfix           integrate\n"
            << "unfix           xlink\n";
    }

    out << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n\n"
        << "write_data      data." << suffix << ".xlink_800 nocoeff\n\n"
        << "# Restore the reacted-bond equilibrium length and 300 K pair interaction\n"
        << "bond_coeff      2 115.4086 " << kBondLengthAngstrom << "\n"
        << "pair_style      lj/gromacs 12 15\n";
    pdms_filler::write_pair_matrix(out, 300.0, false);
    out << "\n";

    if (film) {
        out << "unfix           zlo_wall\n"
            << "unfix           zhi_wall\n"
            << "fix             zlo_wall all wall/lj126 zlo EDGE "
            << cold.epsilon << ' ' << cold.sigma << ' ' << cold.cutoff << " units box\n"
            << "fix             zhi_wall all wall/lj126 zhi EDGE "
            << cold.epsilon << ' ' << cold.sigma << ' ' << cold.cutoff << " units box\n";
    }

    out << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n\n"
        << "# Cool from 800 K to 300 K\n";

    if (film) {
        out << "fix             integrate all npt temp 800.0 300.0 50.0"
            << " x 1.0 1.0 500.0 y 1.0 1.0 500.0 couple xy\n";
    } else {
        out << "fix             integrate all npt temp 800.0 300.0 50.0"
            << " iso 1.0 1.0 500.0\n";
    }

    out << "run             1000000\n"
        << "write_data      data." << suffix << ".300 nocoeff\n"
        << "unfix           integrate\n\n"
        << "# Final 300 K equilibration\n";

    if (film) {
        out << "fix             integrate all npt temp 300.0 300.0 50.0"
            << " x 1.0 1.0 500.0 y 1.0 1.0 500.0 couple xy\n";
    } else {
        out << "fix             integrate all npt temp 300.0 300.0 50.0"
            << " iso 1.0 1.0 500.0\n";
    }

    out << "run             1000000\n"
        << "write_data      data." << suffix << ".npt_eq nocoeff\n\n"
        << "# 1M-step, 300 K NVT trajectory for MSD analysis\n"
        << "unfix           integrate\n"
        << "reset_timestep  0\n\n"
        << "dump            msd all custom " << kMsdDumpEverySteps
        << " dump.msd.lammpstrj id mol type x y z ix iy iz\n"
        << "dump_modify     msd first no sort id\n\n"
        << "fix             msd_integrate all nvt temp 300.0 300.0 50.0\n"
        << "run             " << kMsdProductionSteps << "\n\n"
        << "unfix           msd_integrate\n"
        << "undump          msd\n";

    if (!out) throw std::runtime_error("Failed while writing LAMMPS input file: " + files.input);
}

void write_submit_script(const OutputFiles& files) {
    std::ofstream out(files.submit);
    if (!out) throw std::runtime_error("Cannot open Slurm submit file: " + files.submit);
    out << "#!/bin/bash\n"
        << "#SBATCH --job-name=" << sanitize_job_name(files.case_name) << "_xlink\n"
        << "#SBATCH --time=48:00:00\n"
        << "#SBATCH --nodes=1\n"
        << "#SBATCH --ntasks-per-node=96\n"
        << "#SBATCH --mem=200G\n"
        << "#SBATCH --partition=nova\n"
        << "#SBATCH --mail-user=siteng@iastate.edu\n"
        << "#SBATCH --mail-type=END,FAIL\n"
        << "#SBATCH --output=slurm-%j.out\n"
        << "#SBATCH --error=slurm-%j.err\n\n"
        << "set -euo pipefail\n"
        << "cd -- \"${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}\"\n\n"
        << "module purge\n"
        << "module load intel/22.3.1\n"
        << "module load mpi/2021.7.1\n"
        << "module load lammps/20230802.2-py310-openmpi4-ezoqd7f\n\n"
        << "export OMP_NUM_THREADS=1\n\n"
        << "INPUT=" << shell_single_quote(files.input_basename) << "\n"
        << "OUTPUT=" << shell_single_quote("out." + files.case_name) << "\n"
        << "srun lmp -in \"$INPUT\" > \"$OUTPUT\"\n";
    if (!out) throw std::runtime_error("Failed while writing Slurm submit file: " + files.submit);
}

void write_info(const Settings& s, const System& sys, const Box& box,
                const OutputFiles& files) {
    std::ofstream out(files.info);
    if (!out) throw std::runtime_error("Cannot open model info file: " + files.info);

    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    const char* names[] = {"strands", "crosslinkers", "moderators", "filler"};
    const long long component_beads[] = {
        1LL*s.n1*s.m1,
        1LL*s.n2*s.m2,
        1LL*s.n3*s.m3,
        filler_beads(s)
    };
    const double component_masses[] = {
        component_beads[0] * s.bead_mass,
        component_beads[1] * s.bead_mass,
        component_beads[2] * s.bead_mass,
        s.m4 * filler_chain_mass(s)
    };
    long long total_beads = 0;
    long long total_molecules = 0;
    double total_mass = 0.0;
    for (int i = 0; i < 4; ++i) {
        total_beads += component_beads[i];
        total_molecules += m[i];
        total_mass += component_masses[i];
    }
    const double volume = box.lx * box.ly * box.lz;
    const double accessible_volume = box.lx * box.ly *
        (s.thickness > 0.0 ? s.thickness : box.lz);
    const double box_density = total_mass / (volume * kAvogadroScale);
    const double realized_filler_wt = total_mass == 0.0
        ? 0.0 : 100.0 * component_masses[3] / total_mass;
    const double compression_scale = s.thickness > 0.0
        ? std::sqrt(s.density / s.target_density)
        : std::cbrt(s.density / s.target_density);
    const double initial_to_curing_scale = s.thickness > 0.0
        ? std::sqrt(s.density / kControlledCuringDensity)
        : std::cbrt(s.density / kControlledCuringDensity);
    const double curing_to_target_scale = s.thickness > 0.0
        ? std::sqrt(kControlledCuringDensity / s.target_density)
        : std::cbrt(kControlledCuringDensity / s.target_density);
    const bool controlled_conversion = conversion_control_enabled(s);
    const long long maximum_new_bonds = stoichiometric_maximum_bonds(s);
    const long long requested_new_bonds = target_new_bonds(s);
    const LjParameters hot = lj_parameters(800.0);
    const LjParameters cold = lj_parameters(300.0);
    const std::set<int> sites = crosslinker_sites(s);
    const std::set<int> strand_reactive_sites = strand_sites(s);
    const TopologyCounts strand_counts = strand_topology_counts(s);

    out << std::fixed << std::setprecision(8)
        << "{\n"
        << "  \"format\": \"pdms-elastomer-model-info\",\n"
        << "  \"format_version\": 3,\n"
        << "  \"model\": \"generic PDMS elastomer\",\n"
        << "  \"case_name\": \"" << json_escape(files.case_name) << "\",\n"
        << "  \"geometry\": \"" << (s.thickness > 0.0 ? "film" : "bulk") << "\",\n"
        << "  \"film_thickness_angstrom\": ";
    if (s.thickness > 0.0) out << s.thickness;
    else out << "null";
    out << ",\n  \"film_thickness_source\": ";
    if (s.thickness > 0.0)
        out << "\"requested nominal wall-free thickness at 300 K\"";
    else
        out << "null";
    out << ",\n  \"film_box_Lz_angstrom\": ";
    if (s.thickness > 0.0) out << box.lz;
    else out << "null";
    out << ",\n  \"film_wall_cutoff_reference_temperature_K\": ";
    if (s.thickness > 0.0) out << 300.0;
    else out << "null";
    out << ",\n  \"film_wall_cutoff_per_side_angstrom\": ";
    if (s.thickness > 0.0) out << cold.cutoff;
    else out << "null";
    out << ",\n  \"film_thickness_definition\": ";
    if (s.thickness > 0.0)
        out << "\"Lz - 2*cold_lj.cutoff\"";
    else
        out << "null";
    out << ",\n"
        << "  \"files\": {\n"
        << "    \"data\": \"" << json_escape(files.data_basename) << "\",\n"
        << "    \"lammps_input\": \"" << json_escape(files.input_basename) << "\",\n"
        << "    \"slurm_submit\": \"" << json_escape(files.submit_basename) << "\",\n"
        << "    \"model_info\": \"" << json_escape(files.info_basename) << "\"\n"
        << "  },\n"
        << "  \"generator_input\": {\n"
        << "    \"config_file\": ";
    if (s.config_file.empty()) out << "null";
    else out << '"' << json_escape(s.config_file) << '"';
    out << ",\n"
        << "    \"precedence\": \"defaults < config file < command line\",\n"
        << "    \"resolved\": {\n"
        << "      \"strand_topology\": \"" << s.strand_topology << "\",\n"
        << "      \"strand_length\": ";
    if (s.strand_topology == "star") out << s.star_arm_length;
    else if (s.strand_topology == "grafted") out << "null";
    else out << s.n1;
    out << ",\n"
        << "      \"strand_beads_per_molecule\": " << s.n1 << ",\n"
        << "      \"strand_count\": " << s.m1 << ",\n"
        << "      \"strand_functionality\": " << s.strand_functionality << ",\n"
        << "      \"strand_reactive_distribution\": \""
        << strand_reactive_distribution_name(s) << "\",\n"
        << "      \"star_arm_count\": ";
    if (s.strand_topology == "star") out << s.star_arm_count;
    else out << "null";
    out << ",\n      \"star_center_count\": ";
    if (s.strand_topology == "star")
        out << star_center_count(s.star_arm_count);
    else out << "null";
    out << ",\n      \"star_arm_length\": ";
    if (s.strand_topology == "star") out << s.star_arm_length;
    else out << "null";
    out << ",\n      \"grafted_label\": ";
    if (s.strand_topology == "grafted")
        out << '"' << grafted_architecture_label(s) << '"';
    else out << "null";
    out << ",\n      \"grafted_backbone_length\": ";
    if (s.strand_topology == "grafted") out << s.graft_backbone_length;
    else out << "null";
    out << ",\n      \"grafted_side_chain_length\": ";
    if (s.strand_topology == "grafted") out << s.graft_side_chain_length;
    else out << "null";
    out << ",\n      \"graft_spacing\": ";
    if (s.strand_topology == "grafted") out << s.graft_spacing;
    else out << "null";
    out << ",\n      \"graft_interval\": ";
    if (s.strand_topology == "grafted") out << s.graft_spacing + 1;
    else out << "null";
    out << ",\n      \"grafted_side_chain_count\": ";
    if (s.strand_topology == "grafted") out << s.graft_side_chain_count;
    else out << "null";
    out << ",\n      \"graft_functional_fraction_requested_percent\": ";
    if (s.strand_topology == "grafted")
        out << s.graft_functional_fraction;
    else out << "null";
    out << ",\n      \"graft_functional_side_chain_count\": ";
    if (s.strand_topology == "grafted") out << s.graft_functional_count;
    else out << "null";
    out << ",\n      \"graft_functional_fraction_realized_percent\": ";
    if (s.strand_topology == "grafted")
        out << 100.0*s.graft_functional_count/s.graft_side_chain_count;
    else out << "null";
    out << ",\n"
        << "      \"crosslinker_length\": " << s.n2 << ",\n"
        << "      \"crosslinker_count\": " << s.m2 << ",\n"
        << "      \"crosslinker_functionality\": " << s.crosslinker_functionality << ",\n"
        << "      \"stoichiometry\": \"" << s.stoichiometry_strand_groups
        << ':' << s.stoichiometry_crosslinker_groups << "\",\n"
        << "      \"moderator_count\": " << s.m3 << ",\n"
        << "      \"filler_length\": " << s.n4 << ",\n"
        << "      \"filler_count\": " << s.m4 << ",\n"
        << "      \"requested_filler_weight_percent\": ";
    if (s.m4 > 0) out << s.filler_weight_percent;
    else out << "null";
    out << "\n"
        << "    }\n"
        << "  },\n"
        << "  \"components\": {\n";
    long long first_molecule = 1;
    for (int i = 0; i < 4; ++i) {
        const long long beads = component_beads[i];
        const double weight_percent = total_mass == 0.0
            ? 0.0 : 100.0 * component_masses[i] / total_mass;
        const double molecule_percent = total_molecules == 0 ? 0.0 : 100.0 * m[i] / total_molecules;
        out << "    \"" << names[i] << "\": {\"component\": " << i + 1
            << ", \"N\": " << n[i] << ", \"M\": " << m[i]
            << ", \"beads\": " << beads << ", \"weight_percent\": " << weight_percent
            << ", \"molecule_percent\": " << molecule_percent
            << ", \"molecule_id_start\": ";
        if (m[i] > 0) out << first_molecule;
        else out << "null";
        out << ", \"molecule_id_end\": ";
        if (m[i] > 0) out << first_molecule + m[i] - 1;
        else out << "null";
        out << "}"
            << (i == 3 ? "\n" : ",\n");
        first_molecule += m[i];
    }
    out << "  },\n"
        << "  \"strand\": {\n"
        << "    \"topology\": \"" << s.strand_topology << "\",\n"
        << "    \"beads_per_molecule\": " << s.n1 << ",\n"
        << "    \"functionality\": " << s.strand_functionality << ",\n"
        << "    \"reactive_distribution\": \""
        << strand_reactive_distribution_name(s) << "\",\n"
        << "    \"star_arm_count\": ";
    if (s.strand_topology == "star") out << s.star_arm_count;
    else out << "null";
    out << ",\n    \"star_center_count\": ";
    if (s.strand_topology == "star")
        out << star_center_count(s.star_arm_count);
    else out << "null";
    out << ",\n    \"star_arm_length\": ";
    if (s.strand_topology == "star") out << s.star_arm_length;
    else out << "null";
    out << ",\n    \"grafted_label\": ";
    if (s.strand_topology == "grafted")
        out << '"' << grafted_architecture_label(s) << '"';
    else out << "null";
    out << ",\n    \"grafted_backbone_length\": ";
    if (s.strand_topology == "grafted") out << s.graft_backbone_length;
    else out << "null";
    out << ",\n    \"grafted_side_chain_length\": ";
    if (s.strand_topology == "grafted") out << s.graft_side_chain_length;
    else out << "null";
    out << ",\n    \"graft_spacing\": ";
    if (s.strand_topology == "grafted") out << s.graft_spacing;
    else out << "null";
    out << ",\n    \"graft_interval\": ";
    if (s.strand_topology == "grafted") out << s.graft_spacing + 1;
    else out << "null";
    out << ",\n    \"side_chain_count\": ";
    if (s.strand_topology == "grafted") out << s.graft_side_chain_count;
    else out << "null";
    out << ",\n    \"functional_fraction_requested_percent\": ";
    if (s.strand_topology == "grafted")
        out << s.graft_functional_fraction;
    else out << "null";
    out << ",\n    \"functional_side_chain_count\": ";
    if (s.strand_topology == "grafted") out << s.graft_functional_count;
    else out << "null";
    out << ",\n    \"functional_fraction_realized_percent\": ";
    if (s.strand_topology == "grafted")
        out << 100.0*s.graft_functional_count/s.graft_side_chain_count;
    else out << "null";
    out << ",\n"
        << "    \"reactive_bead_sites\": [";
    bool first_strand_site = true;
    for (int site : strand_reactive_sites) {
        if (!first_strand_site) out << ", ";
        out << site;
        first_strand_site = false;
    }
    out << "],\n"
        << "    \"bonds_per_molecule\": " << strand_counts.bonds << ",\n"
        << "    \"angles_per_molecule\": " << strand_counts.angles << ",\n"
        << "    \"dihedrals_per_molecule\": "
        << strand_counts.dihedrals << "\n"
        << "  },\n"
        << "  \"composition\": {\n"
        << "    \"total_beads\": " << total_beads << ",\n"
        << "    \"hard_maximum_beads\": null,\n"
        << "    \"recommended_maximum_beads\": "
        << kRecommendedTotalBeads << ",\n"
        << "    \"exceeds_recommended_maximum\": "
        << (total_beads > kRecommendedTotalBeads ? "true" : "false")
        << ",\n"
        << "    \"total_mass_g_per_mol_equivalent\": " << total_mass << ",\n"
        << "    \"total_molecules\": " << total_molecules << ",\n"
        << "    \"requested_filler_weight_percent\": ";
    if (s.m4 > 0) out << s.filler_weight_percent;
    else out << "null";
    out << ",\n"
        << "    \"realized_filler_weight_percent\": " << realized_filler_wt << "\n"
        << "  },\n"
        << "  \"pdms_filler\": {\n"
        << "    \"enabled\": " << (s.m4 > 0 ? "true" : "false") << ",\n"
        << "    \"repeat_units_per_chain\": " << s.n4 << ",\n"
        << "    \"chain_count\": " << s.m4 << ",\n"
        << "    \"beads_per_chain\": " << s.n4 << ",\n"
        << "    \"chain_mass_g_per_mol\": " << filler_chain_mass(s) << ",\n"
        << "    \"minimum_separation_angstrom\": " << s.filler_minimum_separation << ",\n"
        << "    \"initial_z_region\": \"central_40_percent_of_box\",\n"
        << "    \"initial_z_lower_fraction_of_Lz\": " << kFillerZLowerFraction << ",\n"
        << "    \"initial_z_upper_fraction_of_Lz\": " << kFillerZUpperFraction << ",\n"
        << "    \"special_bonds_lj\": [0.0, 0.0, 0.5]\n"
        << "  },\n"
        << "  \"force_field\": {\n"
        << "    \"atom_types\": {\n"
        << "      \"1\": {\"name\": \"neutral DMS\", \"mass\": " << s.bead_mass << "},\n"
        << "      \"2\": {\"name\": \"reactive DMS strand/moderator\", \"mass\": " << s.bead_mass << "},\n"
        << "      \"3\": {\"name\": \"reactive DMS crosslinker\", \"mass\": " << s.bead_mass << "}\n"
        << "    },\n"
        << "    \"special_bonds_lj\": [0.0, 0.0, 0.5],\n"
        << "    \"bond_type_map\": {\"ordinary_pdms\": 1, \"moderator_internal\": 1, \"crosslink\": 2}\n"
        << "  },\n"
        << "  \"crosslinker\": {\n"
        << "    \"functionality\": " << s.crosslinker_functionality << ",\n"
        << "    \"distribution\": \"" << json_escape(s.crosslink_distribution) << "\",\n"
        << "    \"reactive_bead_sites\": [";
    bool first = true;
    for (int site : sites) {
        if (!first) out << ", ";
        out << site;
        first = false;
    }
    out << "],\n"
        << "    \"type2_sites_total_including_moderators\": "
        << 1LL*s.strand_functionality*s.m1 + 4LL*s.m3 << ",\n"
        << "    \"type3_sites_total\": " << 1LL*s.m2*s.crosslinker_functionality << "\n"
        << "  },\n"
        << "  \"stoichiometry\": {\n"
        << "    \"definition\": \"strand functional groups : crosslinker functional groups\",\n"
        << "    \"requested\": \"" << s.stoichiometry_strand_groups << ':'
        << s.stoichiometry_crosslinker_groups << "\",\n"
        << "    \"strand_functional_groups\": "
        << 1LL*s.strand_functionality*s.m1 << ",\n"
        << "    \"crosslinker_functional_groups\": "
        << 1LL*s.m2*s.crosslinker_functionality << ",\n"
        << "    \"moderators_included\": false,\n"
        << "    \"extra_moderator_functional_groups\": " << 4LL*s.m3 << "\n"
        << "  },\n"
        << "  \"initial_state\": {\n"
        << "    \"bead_mass_g_per_mol\": " << s.bead_mass << ",\n"
        << "    \"density_g_per_cm3\": " << s.density << ",\n"
        << "    \"density_definition\": \""
        << (s.thickness > 0.0 ?
            "mass / nominal wall-free material volume" : "mass / box volume")
        << "\",\n"
        << "    \"box_density_g_per_cm3\": " << box_density << ",\n"
        << "    \"box_angstrom\": {\"Lx\": " << box.lx << ", \"Ly\": " << box.ly
        << ", \"Lz\": " << box.lz << "},\n"
        << "    \"volume_angstrom3\": " << volume << ",\n"
        << "    \"nominal_material_volume_angstrom3\": "
        << accessible_volume << ",\n"
        << "    \"bond_length_angstrom\": " << kBondLengthAngstrom << ",\n"
        << "    \"placement_spacing_800K_angstrom\": "
        << kPlacementSpacing800KAngstrom << ",\n"
        << "    \"geometry_parameters_source\": \"fixed by the PDMS model\",\n"
        << "    \"intercomponent_minimum_separation_angstrom\": "
        << s.filler_minimum_separation << ",\n"
        << "    \"component_z_placement\": {\n"
        << "      \"strands\": \"bottom_up\",\n"
        << "      \"crosslinker\": \"top_down\",\n"
        << "      \"pdms_filler_fraction_of_Lz\": ["
        << kFillerZLowerFraction << ", " << kFillerZUpperFraction << "],\n"
        << "      \"moderator_fraction_of_Lz\": ["
        << kModeratorZLowerFraction << ", " << kModeratorZUpperFraction << "]\n"
        << "    },\n"
        << "    \"strand_initial_shape\": \""
        << (s.strand_topology == "ring" ? "planar_regular_ring" :
            s.strand_topology == "star" ? "compact_branched_coil" :
            s.strand_topology == "grafted" ?
                "compact_helical_backbone_with_side_chains" :
                                           "straight_linear") << "\"\n"
        << "  },\n"
        << "  \"topology\": {\n"
        << "    \"atoms\": " << sys.atoms.size() << ",\n"
        << "    \"bonds\": " << sys.bonds.size() << ",\n"
        << "    \"angles\": " << sys.angles.size() << ",\n"
        << "    \"dihedrals\": " << sys.dihedrals.size() << "\n"
        << "  },\n"
        << "  \"random_seeds\": {\n"
        << "    \"strand_reactive_site_seed\": " << s.strand_reactive_seed << ",\n"
        << "    \"star_strand_conformation_seed\": " << s.seed << ",\n"
        << "    \"grafted_strand_conformation_seed\": " << s.seed << ",\n"
        << "    \"initial_velocity_seed\": " << s.seed << ",\n"
        << "    \"crosslink_site_seed\": " << s.crosslink_seed << ",\n"
        << "    \"star_moderator_seed\": " << s.seed << ",\n"
        << "    \"pdms_filler_seed\": " << s.filler_seed << ",\n"
        << "    \"bond_creation_seed\": 348154\n"
        << "  },\n"
        << "  \"simulation_template\": {\n";
    if (controlled_conversion) {
        out << "    \"conversion_control\": {\n"
            << "      \"enabled\": true,\n"
            << "      \"requested_percent\": "
            << s.target_conversion_percent << ",\n"
            << "      \"basis\": \"minimum of strand and crosslinker functional groups; moderators excluded\",\n"
            << "      \"strand_functional_groups\": "
            << strand_functional_groups(s) << ",\n"
            << "      \"crosslinker_functional_groups\": "
            << crosslinker_functional_groups(s) << ",\n"
            << "      \"stoichiometric_maximum_new_bonds\": "
            << maximum_new_bonds << ",\n"
            << "      \"target_new_bonds\": " << requested_new_bonds << ",\n"
            << "      \"integer_rounding\": \"floor\",\n"
            << "      \"controller\": \"fix bond/create cumulative count plus fix halt\",\n"
            << "      \"crosslinking_active_during_compression\": false,\n"
            << "      \"curing_density_g_per_cm3\": "
            << kControlledCuringDensity << ",\n"
            << "      \"initial_to_curing_compression_steps\": "
            << kControlledCompressionStageSteps << ",\n"
            << "      \"initial_to_curing_scale_per_deformed_dimension\": "
            << initial_to_curing_scale << ",\n"
            << "      \"conversion_halt_stage\": \"fixed-density 800 K hold at 0.5 g/cm^3\",\n"
            << "      \"halt_check_every_steps\": 1,\n"
            << "      \"maximum_curing_hold_steps\": "
            << kConversionControlledHoldMaximumSteps << ",\n"
            << "      \"curing_to_target_compression_steps\": "
            << kControlledCompressionStageSteps << ",\n"
            << "      \"curing_to_target_scale_per_deformed_dimension\": "
            << curing_to_target_scale << ",\n"
            << "      \"post_cure_800K_equilibration_steps\": "
            << kControlledPostCureEquilibrationSteps << ",\n"
            << "      \"possible_final_step_overshoot\": true\n"
            << "    },\n";
    }
    out << "    \"initial_minimization\": {\n"
        << "      \"enabled\": true,\n"
        << "      \"style\": \"cg\",\n"
        << "      \"energy_tolerance\": 1.0e-4,\n"
        << "      \"force_tolerance_kcal_per_mol_angstrom\": 1.0e-4,\n"
        << "      \"maximum_iterations\": 1000,\n"
        << "      \"maximum_force_evaluations\": 10000,\n"
        << "      \"maximum_atom_displacement_per_iteration_angstrom\": 0.1,\n"
        << "      \"line_search\": \"backtrack\",\n"
        << "      \"film_wall_energy_included\": "
        << (s.thickness > 0.0 ? "true" : "false") << "\n"
        << "    },\n"
        << "    \"target_compressed_density_g_per_cm3\": " << s.target_density << ",\n"
        << "    \"target_density_definition\": \""
        << (s.thickness > 0.0 ?
            "mass / nominal wall-free material volume" : "mass / box volume")
        << "\",\n"
        << "    \"target_box_density_after_deform_g_per_cm3\": "
        << (s.thickness > 0.0 ?
            s.target_density * s.thickness / box.lz : s.target_density)
        << ",\n"
        << "    \"compression_scale_per_deformed_dimension\": " << compression_scale << ",\n"
        << "    \"hot_temperature_K\": 800.0,\n"
        << "    \"initial_velocity_temperature_K\": 800.0,\n"
        << "    \"initial_velocity_distribution\": \"gaussian\",\n"
        << "    \"final_temperature_K\": 300.0,\n"
        << "    \"wall_style\": ";
    if (s.thickness > 0.0)
        out << "\"separate wall/lj126 fixes at zlo EDGE and zhi EDGE\"";
    else out << "null";
    out << ",\n"
        << "    \"hot_lj\": {\"epsilon\": " << hot.epsilon << ", \"sigma\": " << hot.sigma
        << ", \"cutoff\": " << hot.cutoff << "},\n"
        << "    \"cold_lj\": {\"epsilon\": " << cold.epsilon << ", \"sigma\": " << cold.sigma
        << ", \"cutoff\": " << cold.cutoff << "},\n"
        << "    \"timestep_fs\": 5.0,\n"
        << "    \"equilibration_run_steps\": "
        << (controlled_conversion ? 11000000LL : 7000000LL) << ",\n"
        << "    \"total_run_steps\": "
        << (controlled_conversion ? 12000000LL : 8000000LL) << ",\n";
    if (controlled_conversion)
        out << "    \"run_steps_are_upper_bounds\": true,\n";
    out << "    \"bond_creation_active_steps\": "
        << (controlled_conversion ?
            kConversionControlledHoldMaximumSteps : 4000000LL)
        << ",\n"
        << "    \"bond_creation_probability\": "
        << (controlled_conversion ?
            kControlledCrosslinkProbability : s.crosslink_probability) << ",\n"
        << "    \"msd_production\": {\n"
        << "      \"ensemble\": \"NVT\",\n"
        << "      \"temperature_K\": 300.0,\n"
        << "      \"steps\": " << kMsdProductionSteps << ",\n"
        << "      \"duration_ns\": 5.0,\n"
        << "      \"dump_every_steps\": " << kMsdDumpEverySteps << ",\n"
        << "      \"expected_frames\": " << kMsdExpectedFrames << ",\n"
        << "      \"trajectory_file\": \"dump.msd.lammpstrj\",\n"
        << "      \"coordinates\": \"wrapped x y z with image flags ix iy iz\"\n"
        << "    }\n"
        << "  }\n"
        << "}\n";
    if (!out) throw std::runtime_error("Failed while writing model info file: " + files.info);
}

} // namespace

int main(int argc, char** argv) {
    try {
        Settings settings = parse_args(argc, argv);
        resolve_strand_architecture(settings);
        apply_crosslinker_stoichiometry(settings);
        resolve_filler_composition(settings);
        apply_geometry_filename(settings);
        validate(settings);
        report_composition(settings);
        const long long base_beads =
            1LL*settings.n1*settings.m1 + 1LL*settings.n2*settings.m2 +
            1LL*settings.n3*settings.m3;
        const double total_mass =
            base_beads * settings.bead_mass +
            settings.m4 * filler_chain_mass(settings);
        const double volume =
            total_mass / (settings.density * kAvogadroScale);
        Box box{};
        if (settings.thickness > 0.0) {
            // The requested film thickness is the nominal 300 K region beyond
            // both repulsive wall cutoffs. Density is therefore based on the
            // wall-free material volume, not on the two added wall layers.
            const double cold_wall_cutoff = lj_parameters(300.0).cutoff;
            box.lz = settings.thickness + 2.0 * cold_wall_cutoff;
            box.lx = box.ly = std::sqrt(volume / settings.thickness);
        } else {
            box.lx = box.ly = box.lz = std::cbrt(volume);
        }
        const System system = build_system(settings, box);
        const OutputFiles files = output_files(settings);
        create_output_directory(files);
        write_data(settings, system, box, files.data);
        write_lammps_input(settings, files);
        write_submit_script(files);
        write_info(settings, system, box, files);
        std::cerr << "Wrote model package:\n"
                  << "  " << files.data << "\n"
                  << "  " << files.input << "\n"
                  << "  " << files.submit << "\n"
                  << "  " << files.info << "\n"
                  << "System: " << system.atoms.size() << " atoms, "
                  << system.bonds.size() << " bonds, " << system.angles.size() << " angles, "
                  << system.dihedrals.size() << " dihedrals; box "
                  << box.lx << " x " << box.ly << " x " << box.lz << " A\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
