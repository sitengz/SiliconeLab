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
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../Formulation/silicone_oil_component.hpp"
#include "../Formulation/config_input.hpp"

#ifndef SILICONE_FORMULATION_NAME
#define SILICONE_FORMULATION_NAME "V35"
#define SILICONE_DEFAULT_N1 384
#define SILICONE_DEFAULT_M1 306
#define SILICONE_DEFAULT_FUNCTIONALITY 4
#define SILICONE_FOLDED_NETWORK 1
#define SILICONE_CROSSLINK_PROBABILITY 0.5
#endif

namespace {

constexpr double kAvogadroScale = 0.602; // Converts (g/mol)/(g/cm^3) to A^3.
constexpr const char* kFormulationName = SILICONE_FORMULATION_NAME;
constexpr double kOilZLowerFraction = -0.20;
constexpr double kOilZUpperFraction = 0.20;
constexpr double kModeratorZLowerFraction = 0.28;
constexpr double kModeratorZUpperFraction = 0.38;
constexpr long long kMsdProductionSteps = 1000000;
constexpr int kMsdDumpEverySteps = 1000;
constexpr long long kMsdExpectedFrames =
    kMsdProductionSteps / kMsdDumpEverySteps + 1;
constexpr long long kSurfaceRelaxSteps = 10000000;
constexpr long long kSurfaceProductionSteps = 10000000;
constexpr int kSurfaceSampleEverySteps = 1000;
constexpr double kSurfacePaddingAngstrom = 50.0;
constexpr double kCrosslinkCreationProbability =
    SILICONE_CROSSLINK_PROBABILITY;

struct Settings {
    int n1 = SILICONE_DEFAULT_N1, m1 = SILICONE_DEFAULT_M1;
    int n2 = 32,  m2 = 0;   // M2 is resolved from stoichiometry.
    int n3 = 0,   m3 = 0;   // Silicone-oil repeat units and chains.
    int n4 = 5,   m4 = 6;   // Five-bead star-like moderators.
    int crosslinker_functionality = SILICONE_DEFAULT_FUNCTIONALITY;
    std::string crosslink_distribution = "random";
    std::uint32_t crosslink_seed = 20260722u;
    double bead_mass = 74.0;
    double density = 0.1;
    double target_density = 0.8;
    double bond_length = 2.801;
    double spacing = 7.0;
    double thickness = -1.0; // Negative selects the original cubic bulk box.
    std::uint32_t seed = 5489u;
    std::string output = std::string("data.") + kFormulationName + "_no_oil";
    bool output_explicit = false;
    std::string oil = "none";
    bool oil_explicit = false;
    bool oil_length_explicit = false;
    bool oil_weight_explicit = false;
    double oil_weight_percent = -1.0;
    double mps_monomer_percent = -1.0;
    double mps_weight_percent = -1.0;
    bool sequence_explicit = false;
    std::string sequence = "random";
    int mps_per_chain = 0;
    int chains_with_extra_mps = 0;
    long long total_mps_repeats = 0;
    std::string mps_distribution = "fixed";
    std::uint32_t oil_seed = 20260727u;
    double oil_minimum_separation = 4.5;
    std::string config_file;
    double surface_padding = kSurfacePaddingAngstrom;
    int surface_relax_steps = static_cast<int>(kSurfaceRelaxSteps);
    int surface_production_steps = static_cast<int>(kSurfaceProductionSteps);
    int surface_sample_every = kSurfaceSampleEverySteps;
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

void print_help(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << kFormulationName << " four-component formulation generator.\n"
        << "The default contains no component-3 silicone oil.\n\n"
        << "  --n1 N --m1 M    bifunctional network strands\n"
        << "  --n2 N           cross-linker length (default: 32); M2 is stoichiometric\n"
        << "  --n4 N --m4 M    star-moderator size/count (N4 must be 5; defaults: 5, 6)\n"
        << "\nExplicit component-3 oil selection:\n"
        << "  --oil TYPE        pdms, pmps, or copolymer; omit for no oil\n"
        << "  --oil-length N    repeat units per oil chain; required with --oil\n"
        << "  --oil-wt X        oil weight percent of complete formulation; required\n"
        << "  --mps-percent X   MPS repeat-unit percentage for copolymer oil\n"
        << "  --mps-wt X        MPS repeat-unit weight percentage for copolymer oil\n"
        << "  --mps-distribution MODE  fixed (legacy) or balanced across chains\n"
        << "  --sequence MODE   random, alternating, or block (default: random)\n"
        << "  --oil-seed N      oil sequence/conformation/packing seed (default: 20260727)\n"
        << "  --oil-min-separation X  minimum oil/moderator-to-other distance in A"
           " (default: 4.5)\n"
        << "\nNetwork and box controls:\n"
        << "  --functionality F cross-linker reactive sites, 3-16\n"
        << "  --crosslink-distribution MODE\n"
        << "                     reactive-site placement: regular or random (default: random)\n"
        << "  --crosslink-seed N random-site seed (default: 20260722)\n"
        << "  --mass X          bead mass in g/mol (default: 74)\n"
        << "  --density X       initial packing density in g/cm^3 (default: 0.1)\n"
        << "  --target-density X density after scripted compression (default: 0.8)\n"
        << "  --bond-length X   initial bond length in angstrom (default: 2.801)\n"
        << "  --spacing X       spacing between placed molecules (default: 7.0)\n"
        << "  --thickness X     fixed film thickness Lz in angstrom; omit for bulk\n"
        << "  --seed N          reproducible random seed (default: 5489)\n"
        << "  --output FILE     override the automatically generated data filename\n"
        << "                    a case folder is created beside this path\n"
        << "  --config FILE     read key = value settings; CLI values override file\n"
        << "                    film cases also get a separate z-wall-guarded surface job\n"
        << "  --surface-padding X        vacuum per z face in A (default: 50)\n"
        << "  --surface-relax-steps N    surface 300 K NVT steps (default: 10000000)\n"
        << "  --surface-production-steps N  pressure-production NVT steps (default: 10000000)\n"
        << "  --surface-sample-every N   pressure sample interval in steps (default: 1000)\n"
        << "  --help             show this help\n";
}

void apply_option(Settings& s, const std::string& option,
                  const std::string& value) {
        if      (option == "--n1") s.n1 = parse_int(value, option);
        else if (option == "--m1") s.m1 = parse_int(value, option);
        else if (option == "--n2") s.n2 = parse_int(value, option);
        else if (option == "--m2")
            throw std::runtime_error("M2 is determined by M2=2*M1/functionality; do not supply --m2");
        else if (option == "--n3" || option == "--m3" ||
                 option == "--filler-length" || option == "--filler-wt")
            throw std::runtime_error(
                option + " is obsolete; use --oil, --oil-length, and --oil-wt");
        else if (option == "--n4") s.n4 = parse_int(value, option);
        else if (option == "--m4") s.m4 = parse_int(value, option);
        else if (option == "--oil") { s.oil = value; s.oil_explicit = true; }
        else if (option == "--oil-length") {
            s.n3 = parse_int(value, option);
            s.oil_length_explicit = true;
        }
        else if (option == "--oil-wt") {
            s.oil_weight_percent = parse_double(value, option);
            s.oil_weight_explicit = true;
        }
        else if (option == "--mps-percent")
            s.mps_monomer_percent = parse_double(value, option);
        else if (option == "--mps-wt")
            s.mps_weight_percent = parse_double(value, option);
        else if (option == "--mps-distribution") s.mps_distribution = value;
        else if (option == "--sequence") {
            s.sequence = value;
            s.sequence_explicit = true;
        }
        else if (option == "--oil-seed")
            s.oil_seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--oil-min-separation")
            s.oil_minimum_separation = parse_double(value, option);
        else if (option == "--functionality") s.crosslinker_functionality = parse_int(value, option);
        else if (option == "--crosslink-distribution") s.crosslink_distribution = value;
        else if (option == "--crosslink-seed") s.crosslink_seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--mass") s.bead_mass = parse_double(value, option);
        else if (option == "--density") s.density = parse_double(value, option);
        else if (option == "--target-density") s.target_density = parse_double(value, option);
        else if (option == "--bond-length") s.bond_length = parse_double(value, option);
        else if (option == "--spacing") s.spacing = parse_double(value, option);
        else if (option == "--thickness") s.thickness = parse_double(value, option);
        else if (option == "--surface-padding") s.surface_padding = parse_double(value, option);
        else if (option == "--surface-relax-steps") s.surface_relax_steps = parse_int(value, option);
        else if (option == "--surface-production-steps") s.surface_production_steps = parse_int(value, option);
        else if (option == "--surface-sample-every") s.surface_sample_every = parse_int(value, option);
        else if (option == "--seed") s.seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--output") { s.output = value; s.output_explicit = true; }
        else throw std::runtime_error("Unknown option: " + option);
}

Settings parse_args(int argc, char** argv) {
    Settings s;
    s.config_file = silicone_config::parse_arguments(
        argc, argv,
        [&](const std::string& option, const std::string& value) {
            apply_option(s, option, value);
        },
        [&] { print_help(argv[0]); });
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

double oil_chain_mass(const Settings& s) {
    if (s.m3 == 0) return silicone_oil::chain_mass(s.n3, s.mps_per_chain, s.bead_mass);
    return ((1LL * s.m3 * s.n3 - s.total_mps_repeats) * s.bead_mass +
            s.total_mps_repeats * silicone_oil::kMpsRepeatMass) / s.m3;
}

long long oil_beads(const Settings& s) {
    return 1LL * s.m3 * s.n3 + s.total_mps_repeats;
}

std::string oil_name(const Settings& s) {
    if (s.oil == "pdms") return "PDMS";
    if (s.oil == "pmps") return "PMPS";
    if (s.oil == "copolymer") return "Copolymer";
    return "no_oil";
}

void resolve_oil_composition(Settings& s) {
    const bool has_oil_control =
        s.oil_length_explicit || s.oil_weight_explicit ||
        s.mps_monomer_percent >= 0.0 || s.mps_weight_percent >= 0.0 ||
        s.sequence_explicit;
    if (!s.oil_explicit) {
        if (has_oil_control)
            throw std::runtime_error(
                "Oil controls require an explicit --oil pdms, pmps, or copolymer");
        s.oil = "none";
        s.n3 = 0;
        s.m3 = 0;
        s.mps_per_chain = 0;
        if (!s.output_explicit)
            s.output = std::string("data.") + kFormulationName + "_no_oil";
        return;
    }
    if (s.oil != "pdms" && s.oil != "pmps" && s.oil != "copolymer")
        throw std::runtime_error("--oil must be pdms, pmps, or copolymer");
    if (s.mps_distribution != "fixed" && s.mps_distribution != "balanced")
        throw std::runtime_error("--mps-distribution must be fixed or balanced");
    if (s.oil != "copolymer" && s.mps_distribution != "fixed")
        throw std::runtime_error("--mps-distribution balanced requires --oil copolymer");
    if (!s.oil_length_explicit || !s.oil_weight_explicit)
        throw std::runtime_error(
            "--oil requires both --oil-length and --oil-wt");
    if (s.n3 <= 0)
        throw std::runtime_error("--oil-length must be positive");
    if (s.oil_weight_percent <= 0.0 || s.oil_weight_percent >= 100.0)
        throw std::runtime_error("--oil-wt must be greater than 0 and less than 100");
    if (s.oil_minimum_separation <= 0.0 ||
        s.oil_minimum_separation >= 15.0)
        throw std::runtime_error("--oil-min-separation must be greater than 0 and less than 15 A");
    if (s.sequence != "random" && s.sequence != "alternating" &&
        s.sequence != "block")
        throw std::runtime_error("--sequence must be random, alternating, or block");
    if (s.mps_monomer_percent >= 0.0 && s.mps_weight_percent >= 0.0)
        throw std::runtime_error("--mps-percent and --mps-wt are mutually exclusive");

    double requested_mps = 0.0;
    if (s.oil == "pdms" || s.oil == "pmps") {
        if (s.mps_monomer_percent >= 0.0 || s.mps_weight_percent >= 0.0 ||
            s.sequence_explicit)
            throw std::runtime_error(
                "MPS composition and sequence options are only valid with --oil copolymer");
        s.mps_per_chain = s.oil == "pdms" ? 0 : s.n3;
        requested_mps = s.mps_per_chain;
    } else {
        if (s.mps_monomer_percent < 0.0 && s.mps_weight_percent < 0.0)
            throw std::runtime_error(
                "--oil copolymer requires --mps-percent or --mps-wt");
        if (s.mps_weight_percent >= 0.0) {
            if (s.mps_weight_percent <= 0.0 || s.mps_weight_percent >= 100.0)
                throw std::runtime_error("--mps-wt must be greater than 0 and less than 100");
            const double fraction = s.mps_weight_percent / 100.0;
            const double denominator =
                silicone_oil::kMpsRepeatMass * (1.0 - fraction) +
                fraction * s.bead_mass;
            requested_mps =
                fraction * s.n3 * s.bead_mass / denominator;
        } else {
            if (s.mps_monomer_percent <= 0.0 ||
                s.mps_monomer_percent >= 100.0)
                throw std::runtime_error(
                    "--mps-percent must be greater than 0 and less than 100");
            requested_mps =
                s.n3 * s.mps_monomer_percent / 100.0;
        }
        if (s.n3 < 2)
            throw std::runtime_error("Copolymer oil requires --oil-length of at least 2");
        if (s.mps_distribution == "fixed")
            s.mps_per_chain = std::max(
                1, std::min(s.n3 - 1, static_cast<int>(std::lround(requested_mps))));
    }

    const long long base_beads =
        1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n4*s.m4;
    const double base_mass = base_beads * s.bead_mass;
    const double fraction = s.oil_weight_percent / 100.0;
    const double desired_oil_mass = base_mass * fraction / (1.0 - fraction);
    const double estimated_mps = s.mps_distribution == "balanced"
        ? requested_mps : static_cast<double>(s.mps_per_chain);
    const double estimated_chain_mass =
        (s.n3 - estimated_mps) * s.bead_mass +
        estimated_mps * silicone_oil::kMpsRepeatMass;
    s.m3 = std::max(
        1, static_cast<int>(std::lround(desired_oil_mass / estimated_chain_mass)));
    if (s.mps_distribution == "balanced") {
        s.total_mps_repeats = std::max(
            0LL, std::min(1LL * s.n3 * s.m3,
                static_cast<long long>(std::llround(requested_mps * s.m3))));
        s.mps_per_chain = static_cast<int>(s.total_mps_repeats / s.m3);
        s.chains_with_extra_mps = static_cast<int>(s.total_mps_repeats % s.m3);
    } else {
        s.total_mps_repeats = 1LL * s.m3 * s.mps_per_chain;
    }

    if (!s.output_explicit) {
        s.output = std::string("data.") + kFormulationName + '_' +
            oil_name(s) + "_N" + std::to_string(s.n3);
        if (s.oil == "copolymer") {
            if (s.mps_distribution == "balanced")
                s.output += "_MPS" + filename_number(
                    s.mps_weight_percent >= 0.0 ? s.mps_weight_percent : s.mps_monomer_percent) +
                    (s.mps_weight_percent >= 0.0 ? "wt" : "mol") + "_balanced_" + s.sequence;
            else
                s.output += "_MPS" + std::to_string(s.mps_per_chain) + "of" +
                            std::to_string(s.n3) + '_' + s.sequence;
        }
        s.output += "_" + filename_number(s.oil_weight_percent) + "wt";
    }
}

void apply_crosslinker_stoichiometry(Settings& s) {
    if (s.crosslinker_functionality < 3)
        throw std::runtime_error("Cross-linker functionality must be at least 3 to form a network");
    const long long reactive_ends = 2LL * s.m1;
    if (reactive_ends % s.crosslinker_functionality != 0) {
        std::ostringstream message;
        message << "Exact stoichiometry is impossible: 2*M1=" << reactive_ends
                << " is not divisible by functionality=" << s.crosslinker_functionality;
        throw std::runtime_error(message.str());
    }
    const long long molecule_count = reactive_ends / s.crosslinker_functionality;
    if (molecule_count > 100000000LL)
        throw std::runtime_error("Stoichiometric M2 is too large");
    s.m2 = static_cast<int>(molecule_count);
}

void apply_geometry_filename(Settings& s) {
    if (s.thickness > 0.0 && !s.output_explicit)
        s.output += "_film_Lz" + filename_number(s.thickness);
}

void report_composition(const Settings& s) {
    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    const long long component_beads[] = {
        1LL*s.n1*s.m1,
        1LL*s.n2*s.m2,
        oil_beads(s),
        1LL*s.n4*s.m4
    };
    const double component_masses[] = {
        component_beads[0] * s.bead_mass,
        component_beads[1] * s.bead_mass,
        s.m3 * oil_chain_mass(s),
        component_beads[3] * s.bead_mass
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
        if (i == 2 && s.oil != "none")
            std::cerr << ", oil=" << s.oil
                      << ", MPS repeats/chain=" << s.mps_per_chain
                      << (s.chains_with_extra_mps ? " or " + std::to_string(s.mps_per_chain + 1) : "")
                      << ", requested oil wt%=" << s.oil_weight_percent;
        std::cerr << '\n';
    }
    std::cerr << "  total beads=" << total_beads
              << ", total mass=" << total_mass << " g/mol-equivalent\n";
}

void validate(const Settings& s) {
    const int values[] = {s.n1, s.n2, s.n3, s.n4, s.m1, s.m2, s.m3, s.m4};
    if (std::any_of(std::begin(values), std::end(values), [](int x) { return x < 0; }))
        throw std::runtime_error("N and M values cannot be negative");
    if (s.n4 != 5 && s.m4 != 0)
        throw std::runtime_error("The implemented moderator is a five-bead star; use --n4 5");
    if (s.crosslinker_functionality < 3 || s.crosslinker_functionality > 16 || s.crosslinker_functionality > s.n2)
        throw std::runtime_error("Cross-linker functionality must be between 3 and min(16, N2)");
    if (s.crosslink_distribution != "regular" && s.crosslink_distribution != "random")
        throw std::runtime_error("--crosslink-distribution must be regular or random");
    if (s.crosslink_distribution == "regular" && s.crosslinker_functionality > (s.n2 + 1) / 2)
        throw std::runtime_error("Regular placement at 1,3,5,... requires functionality <= ceil(N2/2)");
    if (s.bead_mass <= 0 || s.density <= 0 || s.target_density <= 0 ||
        s.bond_length <= 0 || s.spacing <= 0)
        throw std::runtime_error("Mass, densities, bond length, and spacing must be positive");
    if (s.thickness == 0.0)
        throw std::runtime_error("--thickness must be positive; omit it for the cubic bulk system");
    if (!std::isfinite(s.surface_padding) || s.surface_padding <= 15.0 ||
        s.surface_relax_steps <= 0 || s.surface_production_steps <= 0 ||
        s.surface_sample_every <= 0)
        throw std::runtime_error("Surface padding must exceed the 15 A pair cutoff and surface step counts must be positive");
    const long long total_beads =
        1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + oil_beads(s) + 1LL*s.n4*s.m4;
    if (total_beads <= 0) throw std::runtime_error("The formulation must contain at least one bead");
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

#if SILICONE_FOLDED_NETWORK
struct PathPoint { double x, y; };

std::vector<PathPoint> folded_strand_path(int bead_count, const Box& box,
                                          double bond_length, double spacing) {
    std::vector<PathPoint> path;
    if (bead_count <= 0) return path;
    const double available_x = box.lx - 2.0 * spacing;
    if (available_x < bond_length)
        throw std::runtime_error("Lx is too small to place a folded network strand");
    const int horizontal_steps_max =
        std::max(1, static_cast<int>(std::floor(available_x / bond_length)));
    const int connector_steps_max =
        std::max(1, static_cast<int>(std::ceil(spacing / bond_length)));

    path.reserve(static_cast<size_t>(bead_count));
    double x = 0.0, y = 0.0;
    int direction = 1;
    int horizontal_steps = 0;
    int connector_steps = 0;
    bool connecting_rows = false;

    for (int bead = 0; bead < bead_count; ++bead) {
        path.push_back({x, y});
        if (bead + 1 == bead_count) break;
        if (!connecting_rows && horizontal_steps < horizontal_steps_max) {
            x += direction * bond_length;
            ++horizontal_steps;
        } else if (!connecting_rows) {
            y += bond_length;
            connecting_rows = true;
            connector_steps = 1;
        } else if (connector_steps < connector_steps_max) {
            y += bond_length;
            ++connector_steps;
        } else {
            direction = -direction;
            x += direction * bond_length;
            horizontal_steps = 1;
            connector_steps = 0;
            connecting_rows = false;
        }
    }
    return path;
}

void add_network_strands(System& sys, const Settings& s, const Box& box) {
    if (s.m1 == 0) return;
    if (s.n1 <= 0) throw std::runtime_error("A nonzero M1 requires N1 > 0");
    if (box.ly < 2.0*s.spacing || box.lz < 2.0*s.spacing)
        throw std::runtime_error("Ly and Lz must each be at least twice the placement spacing");

    const std::vector<PathPoint> path =
        folded_strand_path(s.n1, box, s.bond_length, s.spacing);
    double footprint_y = 0.0;
    for (const PathPoint& point : path) footprint_y = std::max(footprint_y, point.y);
    if (footprint_y > box.ly - 2.0*s.spacing)
        throw std::runtime_error("Ly is too small for the folded network-strand footprint");

    const int ny = std::max(1, static_cast<int>(
        std::floor((box.ly - 2.0*s.spacing) / (footprint_y + s.spacing))) + 1);
    const int nz = std::max(1, static_cast<int>(box.lz / s.spacing));
    const long long capacity = 1LL * ny * nz;
    if (capacity < s.m1) {
        std::ostringstream msg;
        msg << "Folded network-strand placement capacity (" << capacity
            << ") is smaller than M1 (" << s.m1 << "). Reduce spacing.";
        throw std::runtime_error(msg.str());
    }

    for (int molecule_index = 0; molecule_index < s.m1; ++molecule_index) {
        const int iy = molecule_index % ny;
        const int iz = molecule_index / ny;
        const bool mirror_x = molecule_index % 2 != 0;
        const double y0 = -box.ly/2 + s.spacing +
                          iy * (footprint_y + s.spacing);
        const double z = -box.lz/2 + s.spacing + iz * s.spacing;
        const int molecule_id =
            static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;

        for (int bead = 1; bead <= s.n1; ++bead) {
            const PathPoint& point = path[static_cast<size_t>(bead - 1)];
            const double x = mirror_x
                ? box.lx/2 - s.spacing - point.x
                : -box.lx/2 + s.spacing + point.x;
            const int atom_type = (bead == 1 || bead == s.n1) ? 2 : 1;
            sys.atoms.push_back({static_cast<int>(sys.atoms.size()) + 1,
                                 molecule_id, atom_type, 0.0, x, y0 + point.y, z});
        }
        add_linear_topology(sys, first_atom, s.n1);
    }
}
#endif

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
    const int nz = component == 3
        ? std::max(1, static_cast<int>((box.lz / 2.0 - spacing) / spacing) + 1)
        : std::max(1, static_cast<int>(box.lz / spacing));
    const long long capacity = 1LL * nx * ny * nz;
    if (capacity < molecule_count) {
        std::ostringstream msg;
        msg << "Placement grid capacity (" << capacity << ") is smaller than M" << component
            << " (" << molecule_count << "). Reduce spacing or molecule length.";
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
        else if (component == 3) z = iz * spacing; // Filler starts at the midplane and grows toward +z.
        else z = -box.lz / 2 + spacing + iz * spacing;

        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 1; bead <= bead_count; ++bead) {
            int atom_type = 1;
            if (component == 1 && (bead == 1 || bead == bead_count)) atom_type = 2;
            if (component == 2 && reactive_sites.count(bead)) atom_type = 3;
            sys.atoms.push_back({static_cast<int>(sys.atoms.size()) + 1, molecule_id, atom_type, 0.0, x, y, z});
            x += reverse ? -bond_length : bond_length;
        }
        add_linear_topology(sys, first_atom, bead_count);
    }
}

double minimum_image(double delta, double length) {
    return delta - std::round(delta / length) * length;
}

bool overlaps_existing(const std::vector<silicone_oil::Vec3>& candidate,
                       const System& sys, const Box& box,
                       double minimum_separation, bool periodic_z) {
    const double minimum_squared = minimum_separation * minimum_separation;
    for (const silicone_oil::Vec3& position : candidate) {
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
    const double margin = s.bond_length;
    if (box.lx <= 2.0*margin || box.ly <= 2.0*margin)
        throw std::runtime_error("Lateral box dimensions are too small for a star moderator");
    std::uniform_real_distribution<double> xdist(-box.lx/2 + margin, box.lx/2 - margin);
    std::uniform_real_distribution<double> ydist(-box.ly/2 + margin, box.ly/2 - margin);
    std::uniform_real_distribution<double> zdist(
        kModeratorZLowerFraction * box.lz,
        kModeratorZUpperFraction * box.lz);
    for (int i = 0; i < s.m4; ++i) {
        std::vector<silicone_oil::Vec3> candidate;
        bool accepted = false;
        for (int attempt = 0; attempt < 10000 && !accepted; ++attempt) {
            const double x = xdist(rng);
            const double y = ydist(rng);
            const double z = zdist(rng);
            const double b = s.bond_length;
            candidate = {
                {x,     y,     z},
                {x + b, y,     z},
                {x - b, y,     z},
                {x,     y + b, z},
                {x,     y - b, z}
            };
            accepted = !overlaps_existing(
                candidate, sys, box, s.oil_minimum_separation,
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
            const silicone_oil::Vec3& position = candidate[bead];
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

void add_silicone_oil(System& sys, const Settings& s, const Box& box) {
    if (s.m3 == 0) return;
    std::vector<silicone_oil::Vec3> existing;
    existing.reserve(sys.atoms.size());
    for (const Atom& atom : sys.atoms)
        existing.push_back({atom.x, atom.y, atom.z});

    silicone_oil::Settings oil;
    oil.length = s.n3;
    oil.chains = s.m3;
    oil.mps_per_chain = s.mps_per_chain;
    oil.chains_with_extra_mps = s.chains_with_extra_mps;
    oil.sequence = s.sequence;
    oil.seed = s.oil_seed;
    oil.minimum_separation = s.oil_minimum_separation;
    oil.z_lower_fraction = kOilZLowerFraction;
    oil.z_upper_fraction = kOilZUpperFraction;
    const silicone_oil::Box oil_box{
        box.lx, box.ly, box.lz, s.thickness <= 0.0
    };
    const silicone_oil::Component generated =
        silicone_oil::generate(oil, oil_box, existing);

    const int atom_offset = static_cast<int>(sys.atoms.size());
    const int molecule_offset =
        sys.atoms.empty() ? 0 : sys.atoms.back().molecule;
    for (const silicone_oil::Atom& atom : generated.atoms) {
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
    for (const silicone_oil::Bond& bond : generated.bonds)
        sys.bonds.push_back({
            static_cast<int>(sys.bonds.size()) + 1, bond.type,
            atom_offset + bond.a, atom_offset + bond.b
        });
    for (const silicone_oil::Angle& angle : generated.angles)
        sys.angles.push_back({
            static_cast<int>(sys.angles.size()) + 1, angle.type,
            atom_offset + angle.a, atom_offset + angle.b,
            atom_offset + angle.c
        });
    for (const silicone_oil::Dihedral& dihedral : generated.dihedrals)
        sys.dihedrals.push_back({
            static_cast<int>(sys.dihedrals.size()) + 1, dihedral.type,
            atom_offset + dihedral.a, atom_offset + dihedral.b,
            atom_offset + dihedral.c, atom_offset + dihedral.d
        });
}

System build_system(const Settings& s, const Box& box) {
    System sys;
    const long long expected_atoms =
        1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + oil_beads(s) + 1LL*s.n4*s.m4;
    if (expected_atoms > 100000000LL) throw std::runtime_error("Requested system is too large");
    sys.atoms.reserve(static_cast<size_t>(expected_atoms));
#if SILICONE_FOLDED_NETWORK
    add_network_strands(sys, s, box);
#else
    add_linear_component(
        sys, s.n1, s.m1, 1, box, s.bond_length, s.spacing);
#endif
    const std::set<int> reactive_sites = crosslinker_sites(s);
    std::cerr << "Cross-linker type-3 sites (" << s.crosslink_distribution << "):";
    for (int site : reactive_sites) std::cerr << ' ' << site;
    std::cerr << '\n';
    add_linear_component(sys, s.n2, s.m2, 2, box, s.bond_length, s.spacing,
                         reactive_sites);
    add_silicone_oil(sys, s, box);
    std::mt19937 rng(s.seed);
    add_star_moderators(sys, s, box, rng);
    if (static_cast<long long>(sys.atoms.size()) != expected_atoms)
        throw std::logic_error("Internal error: generated atom count does not match requested composition");
    return sys;
}

void write_data(const Settings& s, const System& sys, const Box& box,
                const std::string& output_path) {
    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Cannot open output file: " + output_path);
    out << "LAMMPS data file for the " << kFormulationName
        << " formulation with explicit silicone-oil namespace\n\n"
        << sys.atoms.size() << " atoms\n5 atom types\n"
        << sys.bonds.size() << " bonds\n4 bond types\n"
        << sys.angles.size() << " angles\n3 angle types\n"
        << sys.dihedrals.size() << " dihedrals\n4 dihedral types\n\n"
        << std::fixed << std::setprecision(6)
        << -box.lx/2 << ' ' << box.lx/2 << " xlo xhi\n"
        << -box.ly/2 << ' ' << box.ly/2 << " ylo yhi\n"
        << -box.lz/2 << ' ' << box.lz/2 << " zlo zhi\n\nMasses\n\n";
    out << "1 " << s.bead_mass << '\n'
        << "2 " << s.bead_mass << '\n'
        << "3 " << s.bead_mass << '\n'
        << "4 " << silicone_oil::kMpsBackboneMass << '\n'
        << "5 " << silicone_oil::kMpsPendantMass << '\n';
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
    std::string surface_input;
    std::string surface_input_basename;
    std::string surface_submit;
    std::string surface_submit_basename;
    std::string dependent_submit;
    std::string dependent_submit_basename;
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
    if (s.thickness > 0.0) {
        files.surface_input_basename = "in." + files.case_name + ".surface";
        files.surface_submit_basename = "submit." + files.case_name + ".surface.sh";
        files.surface_input = (directory / files.surface_input_basename).string();
        files.surface_submit = (directory / files.surface_submit_basename).string();
        files.dependent_submit_basename = "submit." + files.case_name + ".chain.sh";
        files.dependent_submit = (directory / files.dependent_submit_basename).string();
    }
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
    return result.empty() ? kFormulationName : result;
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
    const LjParameters hot = lj_parameters(800.0);
    const LjParameters cold = lj_parameters(300.0);
    const double hot_global_cutoff =
        silicone_oil::maximum_repulsive_cutoff(800.0);
    const std::string suffix = files.case_name;

    out << "# Generated by the " << kFormulationName << " model generator\n"
        << "# Geometry: " << (film ? "film with fixed Lz" : "periodic bulk") << "\n\n"
        << "units           real\n"
        << "boundary        p p " << (film ? "f" : "p") << "\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     hybrid harmonic quartic\n"
        << "dihedral_style  nharmonic\n"
        << "special_bonds   lj 0 0 0.5\n"
        << std::fixed << std::setprecision(8)
        << "pair_style      lj/cut " << hot_global_cutoff << "\n"
        << "comm_modify     cutoff 15\n"
        << "read_data       " << files.data_basename
        << " extra/bond/per/atom 4 extra/angle/per/atom 10"
        << " extra/dihedral/per/atom 10 extra/special/per/atom 30\n\n"
        << "mass            1 " << s.bead_mass << "\n"
        << "mass            2 " << s.bead_mass << "\n"
        << "mass            3 " << s.bead_mass << "\n"
        << "mass            4 " << silicone_oil::kMpsBackboneMass << "\n"
        << "mass            5 " << silicone_oil::kMpsPendantMass << "\n\n"
        << "bond_coeff      1 115.4086 2.801\n"
        << "bond_coeff      2 115.4086 7.235\n"
        << "bond_coeff      3 108.3835 2.8039\n"
        << "bond_coeff      4 232.8302 3.12497\n"
        << "angle_coeff     1 harmonic 64.62431 111.623\n"
        << "angle_coeff     2 quartic 110.566 64.3974 -139.5241 80.974\n"
        << "angle_coeff     3 quartic 110.746 -20.8906 23.3707 180.3228\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n"
        << "dihedral_coeff  2 8 1.3730 0.2686 0.4017 -1.7250 -0.7052 4.1390 0.2327 -2.2635\n"
        << "dihedral_coeff  3 8 2.3494 -1.5840 -1.6463 3.2133 4.1479 -2.8795 -2.2665 1.1811\n"
        << "dihedral_coeff  4 8 2.23125 0.24735 2.4327 -2.8832 -4.7124 7.17825 2.33495 -3.74\n\n"
        << "# Explicit 800 K repulsive pair matrix\n";
    silicone_oil::write_pair_matrix(out, 800.0, true);
    out << "\n"
        << "neighbor        2 bin\n"
        << "neigh_modify    delay 5 every 1\n"
        << "timestep        5\n"
        << "thermo          1000\n"
        << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n"
        << "restart         100000 restart." << suffix << ".1 restart." << suffix << ".2\n"
        << "dump            traj all custom 100000 dump." << suffix
        << ".lammpstrj id mol type q x y z ix iy iz\n"
        << "dump_modify     traj format line \"%d %d %d %.1f %.3f %.3f %.3f %d %d %d\" sort id\n\n";

    if (film) {
        out << "# Separate repulsive fixes for the lower and upper z box edges.\n"
            << "fix             zlo_wall all wall/lj126 zlo EDGE "
            << hot.epsilon << ' ' << hot.sigma << ' ' << hot.cutoff << " units box\n"
            << "fix             zhi_wall all wall/lj126 zhi EDGE "
            << hot.epsilon << ' ' << hot.sigma << ' ' << hot.cutoff << " units box\n\n";
    }

    out << "# 800 K relaxation\n"
        << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
        << "run             1000000\n"
        << "write_data      data." << suffix << ".rep_800 nocoeff\n\n"
        << "unfix           integrate\n\n";

    if (film) {
        out << "# Lateral compression with crosslinking at fixed film thickness\n"
            << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
            << "fix             xlink all bond/create 1 2 3 " << hot.cutoff
            << " 2 iparam 1 1 jparam 1 1 prob "
            << kCrosslinkCreationProbability << " 348154\n"
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
            << kCrosslinkCreationProbability << " 348154\n"
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
        << "bond_coeff      2 115.4086 2.801\n"
        << "pair_style      lj/gromacs 12 15\n";
    silicone_oil::write_pair_matrix(out, 300.0, false);
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
        << "undump          traj\n"
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

void write_surface_input(const Settings& s, const OutputFiles& files) {
    if (files.surface_input.empty()) return;
    std::ofstream out(files.surface_input);
    if (!out) throw std::runtime_error("Cannot open surface input: " + files.surface_input);
    const std::string& suffix = files.case_name;
    const LjParameters cold = lj_parameters(300.0);
    out << std::fixed << std::setprecision(9)
        << "# Coating surface measurement with distant z guard walls; run after the film's main job.\n"
        << "# The cured film is NOT made by cutting a periodic bulk network.\n"
        << "units           real\n"
        << "boundary        p p f\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     hybrid harmonic quartic\n"
        << "dihedral_style  nharmonic\n"
        << "special_bonds   lj 0 0 0.5\n"
        << "pair_style      lj/gromacs 12 15\n"
        << "comm_modify     cutoff 15\n"
        << "read_data       data." << suffix << ".npt_eq\n\n"
        << "bond_coeff      1 115.4086 2.801\n"
        << "bond_coeff      2 115.4086 2.801\n"
        << "bond_coeff      3 108.3835 2.8039\n"
        << "bond_coeff      4 232.8302 3.12497\n"
        << "angle_coeff     1 harmonic 64.62431 111.623\n"
        << "angle_coeff     2 quartic 110.566 64.3974 -139.5241 80.974\n"
        << "angle_coeff     3 quartic 110.746 -20.8906 23.3707 180.3228\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n"
        << "dihedral_coeff  2 8 1.3730 0.2686 0.4017 -1.7250 -0.7052 4.1390 0.2327 -2.2635\n"
        << "dihedral_coeff  3 8 2.3494 -1.5840 -1.6463 3.2133 4.1479 -2.8795 -2.2665 1.1811\n"
        << "dihedral_coeff  4 8 2.23125 0.24735 2.4327 -2.8832 -4.7124 7.17825 2.33495 -3.74\n";
    silicone_oil::write_pair_matrix(out, 300.0, false);
    out << "\n"
        << "neighbor        2 bin\n"
        << "neigh_modify    delay 5 every 1\n"
        << "timestep        5\n"
        << "# Expand into vacuum without remapping or scaling atom positions.\n"
        << "change_box      all z delta -" << s.surface_padding << ' '
        << s.surface_padding << " units box\n"
        << "# Keep atoms inside the expanded box with separate lower and upper LJ walls.\n"
        << "fix             zlo_wall all wall/lj126 zlo EDGE "
        << cold.epsilon << ' ' << cold.sigma << ' ' << cold.cutoff << " units box\n"
        << "fix             zhi_wall all wall/lj126 zhi EDGE "
        << cold.epsilon << ' ' << cold.sigma << ' ' << cold.cutoff << " units box\n"
        << "thermo          10000\n"
        << "thermo_style    custom step time temp pe pxx pyy pzz lx ly lz\n"
        << "thermo_modify   lost error format float %.12g\n"
        << "variable        surface_time equal time\n"
        << "variable        surface_temp equal temp\n"
        << "variable        surface_pe equal pe\n"
        << "variable        surface_pxx equal pxx\n"
        << "variable        surface_pyy equal pyy\n"
        << "variable        surface_pzz equal pzz\n"
        << "variable        surface_lx equal lx\n"
        << "variable        surface_ly equal ly\n"
        << "variable        surface_lz equal lz\n"
        << "reset_timestep  0\n"
        << "fix             integrate all nvt temp 300.0 300.0 50.0\n"
        << "fix             surface_eq all print " << s.surface_sample_every
        << " \"${surface_time} ${surface_temp} ${surface_pe} ${surface_pxx} ${surface_pyy} "
           "${surface_pzz} ${surface_lx} ${surface_ly} ${surface_lz}\" file stress."
        << suffix << ".surface_eq.dat screen no title \"# time_fs temp_K pe_kcal_per_mol "
           "pxx_atm pyy_atm pzz_atm lx_A ly_A lz_A\"\n"
        << "run             " << s.surface_relax_steps << "\n"
        << "unfix           surface_eq\n"
        << "unfix           integrate\n"
        << "write_data      data." << suffix << ".surface_eq nocoeff\n"
        << "reset_timestep  0\n"
        << "fix             integrate all nvt temp 300.0 300.0 50.0\n"
        << "fix             surface_prod all print " << s.surface_sample_every
        << " \"${surface_time} ${surface_temp} ${surface_pe} ${surface_pxx} ${surface_pyy} "
           "${surface_pzz} ${surface_lx} ${surface_ly} ${surface_lz}\" file stress."
        << suffix << ".surface.dat screen no title \"# time_fs temp_K pe_kcal_per_mol "
           "pxx_atm pyy_atm pzz_atm lx_A ly_A lz_A\"\n"
        << "run             " << s.surface_production_steps << "\n"
        << "unfix           surface_prod\n"
        << "unfix           integrate\n"
        << "write_data      data." << suffix << ".surface_final nocoeff\n";
    if (!out) throw std::runtime_error("Failed while writing surface input: " + files.surface_input);
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

void write_surface_submit_script(const OutputFiles& files) {
    if (files.surface_submit.empty()) return;
    std::ofstream out(files.surface_submit);
    if (!out) throw std::runtime_error("Cannot open surface submit file: " + files.surface_submit);
    out << "#!/bin/bash\n"
        << "#SBATCH --job-name=" << sanitize_job_name(files.case_name) << "_surface\n"
        << "#SBATCH --time=48:00:00\n"
        << "#SBATCH --nodes=1\n"
        << "#SBATCH --ntasks-per-node=96\n"
        << "#SBATCH --mem=200G\n"
        << "#SBATCH --partition=nova\n"
        << "#SBATCH --mail-user=siteng@iastate.edu\n"
        << "#SBATCH --mail-type=END,FAIL\n"
        << "#SBATCH --output=slurm-%j.surface.out\n"
        << "#SBATCH --error=slurm-%j.surface.err\n\n"
        << "set -euo pipefail\n"
        << "cd -- \"${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}\"\n"
        << "test -s " << shell_single_quote("data." + files.case_name + ".npt_eq")
        << " || { echo 'Run the film curing job first' >&2; exit 1; }\n"
        << "module purge\n"
        << "module load intel/22.3.1\n"
        << "module load mpi/2021.7.1\n"
        << "module load lammps/20230802.2-py310-openmpi4-ezoqd7f\n"
        << "export OMP_NUM_THREADS=1\n"
        << "INPUT=" << shell_single_quote(files.surface_input_basename) << "\n"
        << "OUTPUT=" << shell_single_quote("out." + files.case_name + ".surface") << "\n"
        << "srun lmp -in \"$INPUT\" > \"$OUTPUT\"\n";
    if (!out) throw std::runtime_error("Failed while writing surface submit file: " + files.surface_submit);
}

void write_dependent_submit_script(const OutputFiles& files) {
    if (files.dependent_submit.empty()) return;
    std::ofstream out(files.dependent_submit);
    if (!out) throw std::runtime_error("Cannot open dependent submit file: " + files.dependent_submit);
    out << "#!/usr/bin/env bash\n"
        << "# Submit the film job now and the surface job after it succeeds.\n"
        << "# This is a login-node launcher, not a Slurm batch job.\n"
        << "set -euo pipefail\n"
        << "cd -- \"$(dirname -- \"$0\")\"\n"
        << "FILM_SCRIPT=" << shell_single_quote(files.submit_basename) << "\n"
        << "SURFACE_SCRIPT=" << shell_single_quote(files.surface_submit_basename) << "\n"
        << "for required in " << shell_single_quote(files.data_basename)
        << " " << shell_single_quote(files.input_basename)
        << " " << shell_single_quote(files.surface_input_basename)
        << " \"$FILM_SCRIPT\" \"$SURFACE_SCRIPT\"; do\n"
        << "    [[ -s $required ]] || { printf 'Missing or empty: %s\\n' \"$required\" >&2; exit 1; }\n"
        << "done\n"
        << "if [[ ${1-} == --dry-run ]]; then\n"
        << "    printf 'Would submit %s, then %s with afterok dependency.\\n' \"$FILM_SCRIPT\" \"$SURFACE_SCRIPT\"\n"
        << "    exit 0\n"
        << "fi\n"
        << "[[ $# -eq 0 ]] || { echo 'Usage: bash submit.<case>.chain.sh [--dry-run]' >&2; exit 2; }\n"
        << "command -v sbatch >/dev/null || { echo 'sbatch is not available' >&2; exit 1; }\n"
        << "film_submission=$(sbatch --parsable \"$FILM_SCRIPT\")\n"
        << "film_id=${film_submission%%;*}\n"
        << "[[ $film_id =~ ^[0-9]+$ ]] || { printf 'Unexpected film job ID: %s\\n' \"$film_submission\" >&2; exit 1; }\n"
        << "printf 'Film job submitted: %s\\n' \"$film_id\"\n"
        << "if ! surface_submission=$(sbatch --parsable --dependency=\"afterok:$film_id\" \"$SURFACE_SCRIPT\"); then\n"
        << "    printf 'Surface submission failed; film job %s is already queued.\\n' \"$film_id\" >&2\n"
        << "    exit 1\n"
        << "fi\n"
        << "surface_id=${surface_submission%%;*}\n"
        << "[[ $surface_id =~ ^[0-9]+$ ]] || { printf 'Unexpected surface job ID: %s\\n' \"$surface_submission\" >&2; exit 1; }\n"
        << "printf 'Surface job submitted: %s (afterok:%s)\\n' \"$surface_id\" \"$film_id\"\n";
    if (!out) throw std::runtime_error("Failed while writing dependent submit file: " + files.dependent_submit);
}
void write_info(const Settings& s, const System& sys, const Box& box,
                const OutputFiles& files) {
    std::ofstream out(files.info);
    if (!out) throw std::runtime_error("Cannot open model info file: " + files.info);

    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    const char* names[] = {"network_strands", "crosslinkers", "silicone_oil_filler",
                           "star_moderators"};
    const long long component_beads[] = {
        1LL*s.n1*s.m1,
        1LL*s.n2*s.m2,
        oil_beads(s),
        1LL*s.n4*s.m4
    };
    const double component_masses[] = {
        component_beads[0] * s.bead_mass,
        component_beads[1] * s.bead_mass,
        s.m3 * oil_chain_mass(s),
        component_beads[3] * s.bead_mass
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
    const double realized_filler_wt = total_mass == 0.0
        ? 0.0 : 100.0 * component_masses[2] / total_mass;
    const double compression_scale = s.thickness > 0.0
        ? std::sqrt(s.density / s.target_density)
        : std::cbrt(s.density / s.target_density);
    const LjParameters hot = lj_parameters(800.0);
    const LjParameters cold = lj_parameters(300.0);
    const std::set<int> sites = crosslinker_sites(s);

    out << std::fixed << std::setprecision(8)
        << "{\n"
        << "  \"format\": \"" << kFormulationName << "-model-info\",\n"
        << "  \"format_version\": 2,\n"
        << "  \"formulation\": \"" << kFormulationName << "\",\n"
        << "  \"case_name\": \"" << json_escape(files.case_name) << "\",\n"
        << "  \"geometry\": \"" << (s.thickness > 0.0 ? "film" : "bulk") << "\",\n"
        << "  \"film_thickness_angstrom\": ";
    if (s.thickness > 0.0) out << s.thickness;
    else out << "null";
    out << ",\n  \"film_thickness_source\": ";
    if (s.thickness > 0.0)
        out << "\"requested value captured from the equilibrated bulk result\"";
    else
        out << "null";
    out << ",\n"
        << "  \"files\": {\n"
        << "    \"data\": \"" << json_escape(files.data_basename) << "\",\n"
        << "    \"lammps_input\": \"" << json_escape(files.input_basename) << "\",\n"
        << "    \"slurm_submit\": \"" << json_escape(files.submit_basename) << "\",\n"
        << "    \"model_info\": \"" << json_escape(files.info_basename) << "\",\n"
        << "    \"surface_input\": ";
    if (files.surface_input_basename.empty()) out << "null";
    else out << '"' << json_escape(files.surface_input_basename) << '"';
    out << ",\n    \"surface_submit\": ";
    if (files.surface_submit_basename.empty()) out << "null";
    else out << '"' << json_escape(files.surface_submit_basename) << '"';
    out << ",\n    \"dependent_submit\": ";
    if (files.dependent_submit_basename.empty()) out << "null";
    else out << '"' << json_escape(files.dependent_submit_basename) << '"';
    out << "\n"
        << "  },\n"
        << "  \"generator_input\": {\n"
        << "    \"config_file\": ";
    if (s.config_file.empty()) out << "null";
    else out << '"' << json_escape(s.config_file) << '"';
    out << ",\n"
        << "    \"precedence\": \"defaults < config file < command line\"\n"
        << "  },\n"
        << "  \"components\": {\n";
    for (int i = 0; i < 4; ++i) {
        const long long beads = component_beads[i];
        const double weight_percent = total_mass == 0.0
            ? 0.0 : 100.0 * component_masses[i] / total_mass;
        const double molecule_percent = total_molecules == 0 ? 0.0 : 100.0 * m[i] / total_molecules;
        out << "    \"" << names[i] << "\": {\"N\": " << n[i] << ", \"M\": " << m[i]
            << ", \"beads\": " << beads << ", \"weight_percent\": " << weight_percent
            << ", \"molecule_percent\": " << molecule_percent << "}"
            << (i == 3 ? "\n" : ",\n");
    }
    out << "  },\n"
        << "  \"composition\": {\n"
        << "    \"total_beads\": " << total_beads << ",\n"
        << "    \"total_mass_g_per_mol_equivalent\": " << total_mass << ",\n"
        << "    \"total_molecules\": " << total_molecules << ",\n"
        << "    \"requested_filler_weight_percent\": ";
    if (s.oil != "none") out << s.oil_weight_percent;
    else out << "null";
    out << ",\n"
        << "    \"realized_filler_weight_percent\": " << realized_filler_wt << "\n"
        << "  },\n"
        << "  \"silicone_oil\": {\n"
        << "    \"model\": \"" << s.oil << "\",\n"
        << "    \"explicitly_selected\": " << (s.oil_explicit ? "true" : "false") << ",\n"
        << "    \"repeat_units_per_chain\": " << s.n3 << ",\n"
        << "    \"chain_count\": " << s.m3 << ",\n"
        << "    \"dms_repeats_per_chain\": " << s.n3 - s.mps_per_chain << ",\n"
        << "    \"mps_repeats_per_chain\": " << s.mps_per_chain << ",\n"
        << "    \"mps_distribution\": \"" << json_escape(s.mps_distribution) << "\",\n"
        << "    \"chains_with_extra_mps\": " << s.chains_with_extra_mps << ",\n"
        << "    \"total_mps_repeats\": " << s.total_mps_repeats << ",\n"
        << "    \"realized_mps_mol_percent\": "
        << (s.m3 == 0 ? 0.0 : 100.0 * s.total_mps_repeats / (1LL * s.m3 * s.n3)) << ",\n"
        << "    \"beads_per_chain\": " << s.n3 + s.mps_per_chain << ",\n"
        << "    \"average_beads_per_chain\": "
        << (s.m3 == 0 ? 0.0 : static_cast<double>(oil_beads(s)) / s.m3) << ",\n"
        << "    \"chain_mass_g_per_mol\": " << oil_chain_mass(s) << ",\n"
        << "    \"sequence\": \"" << json_escape(s.sequence) << "\",\n"
        << "    \"minimum_separation_angstrom\": " << s.oil_minimum_separation << ",\n"
        << "    \"initial_z_region\": \"central_40_percent_of_box\",\n"
        << "    \"initial_z_lower_fraction_of_Lz\": " << kOilZLowerFraction << ",\n"
        << "    \"initial_z_upper_fraction_of_Lz\": " << kOilZUpperFraction << ",\n"
        << "    \"special_bonds_lj\": [0.0, 0.0, 0.5]\n"
        << "  },\n"
        << "  \"force_field\": {\n"
        << "    \"atom_types\": {\n"
        << "      \"1\": {\"name\": \"neutral DMS\", \"mass\": " << s.bead_mass << "},\n"
        << "      \"2\": {\"name\": \"reactive DMS strand/moderator\", \"mass\": " << s.bead_mass << "},\n"
        << "      \"3\": {\"name\": \"reactive DMS crosslinker\", \"mass\": " << s.bead_mass << "},\n"
        << "      \"4\": {\"name\": \"MPS backbone\", \"mass\": " << silicone_oil::kMpsBackboneMass << "},\n"
        << "      \"5\": {\"name\": \"MPS pendant\", \"mass\": " << silicone_oil::kMpsPendantMass << "}\n"
        << "    },\n"
        << "    \"special_bonds_lj\": [0.0, 0.0, 0.5],\n"
        << "    \"dms_mps_mixing\": \"0.579966 * geometric epsilon and arithmetic sigma\",\n"
        << "    \"dms_mps_epsilon_factor\": "
        << silicone_oil::kDmsMpsEpsilonFactor << ",\n"
        << "    \"oil_bond_type_remap\": {\"dms_or_mixed\": 1, \"moderator_internal\": 1, \"formulation_crosslink\": 2, \"mps_backbone\": 3, \"mps_pendant\": 4}\n"
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
        << "    \"type2_sites_total_including_moderators\": " << 2LL*s.m1 + 4LL*s.m4 << ",\n"
        << "    \"type3_sites_total\": " << 1LL*s.m2*s.crosslinker_functionality << "\n"
        << "  },\n"
        << "  \"initial_state\": {\n"
        << "    \"bead_mass_g_per_mol\": " << s.bead_mass << ",\n"
        << "    \"density_g_per_cm3\": " << s.density << ",\n"
        << "    \"box_angstrom\": {\"Lx\": " << box.lx << ", \"Ly\": " << box.ly
        << ", \"Lz\": " << box.lz << "},\n"
        << "    \"volume_angstrom3\": " << volume << ",\n"
        << "    \"bond_length_angstrom\": " << s.bond_length << ",\n"
        << "    \"placement_spacing_angstrom\": " << s.spacing << ",\n"
        << "    \"intercomponent_minimum_separation_angstrom\": "
        << s.oil_minimum_separation << ",\n"
        << "    \"component_z_placement\": {\n"
        << "      \"network\": \"bottom_up\",\n"
        << "      \"crosslinker\": \"top_down\",\n"
        << "      \"silicone_oil_fraction_of_Lz\": ["
        << kOilZLowerFraction << ", " << kOilZUpperFraction << "],\n"
        << "      \"moderator_fraction_of_Lz\": ["
        << kModeratorZLowerFraction << ", " << kModeratorZUpperFraction << "]\n"
        << "    },\n"
        << "    \"network_strand_initial_shape\": \""
        << (SILICONE_FOLDED_NETWORK ? "folded_serpentine" : "straight_linear")
        << "\"\n"
        << "  },\n"
        << "  \"topology\": {\n"
        << "    \"atoms\": " << sys.atoms.size() << ",\n"
        << "    \"bonds\": " << sys.bonds.size() << ",\n"
        << "    \"angles\": " << sys.angles.size() << ",\n"
        << "    \"dihedrals\": " << sys.dihedrals.size() << "\n"
        << "  },\n"
        << "  \"random_seeds\": {\n"
        << "    \"crosslink_site_seed\": " << s.crosslink_seed << ",\n"
        << "    \"star_moderator_seed\": " << s.seed << ",\n"
        << "    \"silicone_oil_seed\": " << s.oil_seed << ",\n"
        << "    \"bond_creation_seed\": 348154\n"
        << "  },\n"
        << "  \"simulation_template\": {\n"
        << "    \"target_compressed_density_g_per_cm3\": " << s.target_density << ",\n"
        << "    \"compression_scale_per_deformed_dimension\": " << compression_scale << ",\n"
        << "    \"hot_temperature_K\": 800.0,\n"
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
        << "    \"equilibration_run_steps\": 7000000,\n"
        << "    \"total_run_steps\": 8000000,\n"
        << "    \"bond_creation_active_steps\": 4000000,\n"
        << "    \"bond_creation_probability\": "
        << kCrosslinkCreationProbability << ",\n"
        << "    \"surface_measurement\": {\n"
        << "      \"enabled\": " << (s.thickness > 0.0 ? "true" : "false") << ",\n"
        << "      \"input_data\": ";
    if (s.thickness > 0.0) out << '"' << "data." << json_escape(files.case_name) << ".npt_eq\"";
    else out << "null";
    out << ",\n"
        << "      \"vacuum_padding_per_face_angstrom\": " << s.surface_padding << ",\n"
        << "      \"z_guard_walls\": {\"style\": \"wall/lj126\", \"locations\": [\"zlo EDGE\", \"zhi EDGE\"], \"epsilon_kcal_per_mol\": " << cold.epsilon << ", \"sigma_angstrom\": " << cold.sigma << ", \"cutoff_angstrom\": " << cold.cutoff << "},\n"
        << "      \"ensemble\": \"NVT\",\n"
        << "      \"temperature_K\": 300.0,\n"
        << "      \"relaxation_steps\": " << s.surface_relax_steps << ",\n"
        << "      \"production_steps\": " << s.surface_production_steps << ",\n"
        << "      \"sample_every_steps\": " << s.surface_sample_every << ",\n"
        << "      \"observable\": \"apparent mechanical surface stress, not surface free energy\"\n"
        << "    },\n"
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
        apply_crosslinker_stoichiometry(settings);
        resolve_oil_composition(settings);
        apply_geometry_filename(settings);
        validate(settings);
        report_composition(settings);
        const long long base_beads =
            1LL*settings.n1*settings.m1 + 1LL*settings.n2*settings.m2 +
            1LL*settings.n4*settings.m4;
        const double total_mass =
            base_beads * settings.bead_mass +
            settings.m3 * oil_chain_mass(settings);
        const double volume =
            total_mass / (settings.density * kAvogadroScale);
        Box box{};
        if (settings.thickness > 0.0) {
            box.lz = settings.thickness;
            box.lx = box.ly = std::sqrt(volume / box.lz);
        } else {
            box.lx = box.ly = box.lz = std::cbrt(volume);
        }
        const System system = build_system(settings, box);
        const OutputFiles files = output_files(settings);
        create_output_directory(files);
        write_data(settings, system, box, files.data);
        write_lammps_input(settings, files);
        write_submit_script(files);
        write_surface_input(settings, files);
        write_surface_submit_script(files);
        write_dependent_submit_script(files);
        write_info(settings, system, box, files);
        std::cerr << "Wrote model package:\n"
                  << "  " << files.data << "\n"
                  << "  " << files.input << "\n"
                  << "  " << files.submit << "\n"
                  << "  " << files.info << "\n"
                  << (files.surface_input.empty() ? "" : "  " + files.surface_input + "\n")
                  << (files.surface_submit.empty() ? "" : "  " + files.surface_submit + "\n")
                  << (files.dependent_submit.empty() ? "" : "  " + files.dependent_submit + "\n")
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
