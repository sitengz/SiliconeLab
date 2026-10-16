#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "config_input.hpp"

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kAvogadroScale = 0.602214076;
constexpr double kBoundaryClearance = 1.0e-4;
constexpr double kTimestepFs = 5.0;
constexpr long long kBulkFinalNptSteps = 5000000;
constexpr long long kEquilibrationSteps = 1000000 + 1000000 + 3000000 + 1000000 + kBulkFinalNptSteps;
constexpr int kEnergySampleEverySteps = 1000;
constexpr long long kFilmWallSteps = 100000;
constexpr long long kFilmRelaxSteps = 10000000;
constexpr long long kFilmProductionSteps = 10000000;
constexpr double kMinimumFilmPaddingAngstrom = 20.0;

constexpr double kDmsMass = 74.0;
constexpr double kMpsBackboneMass = 59.1204;
constexpr double kMpsPendantMass = 77.106;
constexpr double kMpsRepeatMass = kMpsBackboneMass + kMpsPendantMass;

constexpr double kDmsBondLength = 2.801;
constexpr double kMpsBackboneBondLength = 2.8039;
constexpr double kMpsPendantBondLength = 3.12497;
constexpr double kDmsAngleDegrees = 111.623;
constexpr double kMpsBackboneAngleDegrees = 110.566;
constexpr double kMpsPendantAngleDegrees = 110.746;

constexpr double kDmsEpsilon = 1.01287845;
constexpr double kDmsSigma = 6.44584366;
constexpr double kMpsBackboneEpsilon = 1.779864125;
constexpr double kMpsBackboneSigma = 5.475173;
constexpr double kMpsPendantEpsilon = 0.848858275;
constexpr double kMpsPendantSigma = 5.96750575;
constexpr double kMpsBackbonePendantEpsilon = 1.33135725;
constexpr double kMpsBackbonePendantSigma = 5.72729065;
constexpr double kDmsMpsEpsilonFactor = 0.579966;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(double scale, const Vec3& v) {
    return {scale * v.x, scale * v.y, scale * v.z};
}

Vec3 operator/(const Vec3& v, double scale) {
    return {v.x / scale, v.y / scale, v.z / scale};
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double norm2(const Vec3& v) {
    return dot(v, v);
}

double norm(const Vec3& v) {
    return std::sqrt(norm2(v));
}

Vec3 normalized(const Vec3& v) {
    const double length = norm(v);
    if (length < 1.0e-12)
        throw std::runtime_error("Cannot normalize a zero-length vector");
    return v / length;
}

double radians(double degrees) {
    return degrees * kPi / 180.0;
}

template <typename T>
T clamp_value(T value, T lower, T upper) {
    return std::max(lower, std::min(value, upper));
}

struct Settings {
    int length = 16;
    int chains = 625;
    bool length_explicit = false;
    bool chains_explicit = false;
    bool length_list_mode = false;
    std::vector<std::pair<int, int>> chain_count_entries;
    std::vector<int> molecule_lengths;
    std::vector<int> molecule_mps_counts;
    double mps_monomer_percent = 100.0;
    double mps_weight_percent = -1.0;
    bool mps_percent_explicit = false;
    std::string sequence = "random";
    double density = 0.1;
    double target_density = 0.8;
    double film_padding = -1.0; // Negative selects the 20 A safety minimum.
    double minimum_separation = 4.5;
    std::uint32_t seed = 20260727u;
    std::uint32_t velocity_seed = 492845u;
    std::string output;
    bool output_explicit = false;
    long long total_repeats = 0;
    long long total_mps_repeats = 0;
    int base_mps_per_chain = 0;
    int chains_with_extra_mps = 0;
    std::string config_file;
};

double film_padding_angstrom(const Settings& settings) {
    return std::max(kMinimumFilmPaddingAngstrom, settings.film_padding);
}

struct Atom {
    int id = 0;
    int molecule = 0;
    int type = 0;
    double charge = 0.0;
    Vec3 position;
    int ix = 0;
    int iy = 0;
    int iz = 0;
};

struct Bond {
    int id, type, a, b;
};

struct Angle {
    int id, type, a, b, c;
};

struct Dihedral {
    int id, type, a, b, c, d;
};

struct System {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
};

struct Box {
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
};

struct PairParameters {
    double epsilon;
    double sigma;
};

double maximum_repulsive_cutoff(double temperature);

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
    std::string film_equilibration_energy_basename;
    std::string film_energy_basename;
    std::string film_input;
    std::string film_input_basename;
    std::string film_submit;
    std::string film_submit_basename;
    std::string pair_submit;
    std::string pair_submit_basename;
    std::string bulk_equilibrated_data_basename;
    std::string film_initial_data_basename;
    std::string film_equilibrated_data_basename;
    std::string film_final_data_basename;
    std::string unwrapped_dump_basename;
    std::string case_name;
};

struct LocalMolecule {
    std::vector<Vec3> backbone;
    std::vector<Vec3> pendants;
};

struct LocalSite {
    Vec3 position;
    int monomer;
    bool pendant;
};

int parse_int(const std::string& value, const std::string& option) {
    try {
        std::size_t used = 0;
        const int result = std::stoi(value, &used);
        if (used != value.size()) throw std::invalid_argument("");
        return result;
    } catch (...) {
        throw std::runtime_error("Invalid integer for " + option + ": " + value);
    }
}

double parse_double(const std::string& value, const std::string& option) {
    try {
        std::size_t used = 0;
        const double result = std::stod(value, &used);
        if (used != value.size() || !std::isfinite(result)) throw std::invalid_argument("");
        return result;
    } catch (...) {
        throw std::runtime_error("Invalid number for " + option + ": " + value);
    }
}

void print_help(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "Standalone PDMS/PMPS silicone-oil generator.\n"
        << "Each run writes bulk data and a dependent film workflow.\n"
        << "One DMS repeat is one type-1 bead. One MPS repeat is a type-4\n"
        << "backbone bead with one type-5 pendant bead.\n\n"
        << "  --length N              repeat units per oil chain (default: 16)\n"
        << "  --chains M              number of oil chains (default: 625)\n"
        << "  --chain-count 'N COUNT'  repeatable length/count row; in a config, use\n"
        << "                           chain_count = N COUNT (cannot mix with length/chains)\n"
        << "  --n N                    alias for --length\n"
        << "  --m M                    alias for --chains\n"
        << "  --mps-percent X          overall MPS monomer percentage, 0-100 (default: 100)\n"
        << "  --mps-wt X               target MPS repeat-unit weight percentage, 0-100\n"
        << "                           (mutually exclusive with --mps-percent)\n"
        << "  --sequence MODE          random, alternating, or block (default: random)\n"
        << "  --density X              initial mass density in g/cm^3 (default: 0.1)\n"
        << "  --target-density X       density after 800 K compression (default: 0.8)\n"
        << "  --film-padding X         requested vacuum at each z face in A\n"
        << "                           (minimum/default: 20 A; smaller values raised to 20 A)\n"
        << "  --min-separation X       minimum intermolecular distance in A (default: 4.5)\n"
        << "  --seed N                 conformation/packing seed (default: 20260727)\n"
        << "  --velocity-seed N        LAMMPS velocity seed (default: 492845)\n"
        << "  --output FILE            name the initial data file; with --config,\n"
        << "                           relative paths resolve beside the config\n"
        << "  --config FILE            read key = value settings; CLI values override file\n"
        << "  --help                   show this help\n";
}

void apply_option(Settings& settings, const std::string& option,
                  const std::string& value) {
    if (option == "--length" || option == "--n") {
        settings.length = parse_int(value, option);
        settings.length_explicit = true;
    } else if (option == "--chains" || option == "--m") {
        settings.chains = parse_int(value, option);
        settings.chains_explicit = true;
    } else if (option == "--chain-count") {
        std::istringstream row(value);
        std::string length_text, count_text, extra;
        if (!(row >> length_text >> count_text) || (row >> extra))
            throw std::runtime_error("--chain-count needs two integers: length count");
        settings.chain_count_entries.emplace_back(
            parse_int(length_text, option), parse_int(count_text, option));
    } else if (option == "--mps-percent") {
        settings.mps_monomer_percent = parse_double(value, option);
        settings.mps_percent_explicit = true;
    } else if (option == "--mps-wt")
        settings.mps_weight_percent = parse_double(value, option);
    else if (option == "--sequence")
        settings.sequence = value;
    else if (option == "--density")
        settings.density = parse_double(value, option);
    else if (option == "--target-density")
        settings.target_density = parse_double(value, option);
    else if (option == "--film-padding")
        settings.film_padding = parse_double(value, option);
    else if (option == "--min-separation")
        settings.minimum_separation = parse_double(value, option);
    else if (option == "--seed") {
        const int parsed = parse_int(value, option);
        if (parsed <= 0) throw std::runtime_error("--seed must be positive");
        settings.seed = static_cast<std::uint32_t>(parsed);
    } else if (option == "--velocity-seed") {
        const int parsed = parse_int(value, option);
        if (parsed <= 0) throw std::runtime_error("--velocity-seed must be positive");
        settings.velocity_seed = static_cast<std::uint32_t>(parsed);
    } else if (option == "--output") {
        settings.output = value;
        settings.output_explicit = true;
    } else {
        throw std::runtime_error("Unknown option: " + option);
    }
}

Settings parse_args(int argc, char** argv) {
    Settings settings;
    settings.config_file = silicone_config::parse_arguments(
        argc, argv,
        [&](const std::string& option, const std::string& value) {
            apply_option(settings, option, value);
        },
        [&] { print_help(argv[0]); });
    return settings;
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

void resolve_composition(Settings& settings) {
    settings.length_list_mode = !settings.chain_count_entries.empty();
    if (settings.length_list_mode) {
        if (settings.length_explicit || settings.chains_explicit)
            throw std::runtime_error("--chain-count cannot be mixed with --length or --chains");
        std::sort(settings.chain_count_entries.begin(), settings.chain_count_entries.end());
        long long chain_total = 0;
        long long repeat_total = 0;
        int previous_length = 0;
        for (const auto& entry : settings.chain_count_entries) {
            const int length = entry.first;
            const int count = entry.second;
            if (length <= 0 || count <= 0)
                throw std::runtime_error("--chain-count lengths and counts must be positive");
            if (length == previous_length)
                throw std::runtime_error("--chain-count contains a duplicate length");
            previous_length = length;
            chain_total += count;
            repeat_total += 1LL * length * count;
            if (chain_total > std::numeric_limits<int>::max() ||
                repeat_total > std::numeric_limits<int>::max())
                throw std::runtime_error("--chain-count exceeds 32-bit LAMMPS limits");
        }
        settings.chains = static_cast<int>(chain_total);
        settings.total_repeats = repeat_total;
        settings.molecule_lengths.reserve(static_cast<std::size_t>(settings.chains));
        for (const auto& entry : settings.chain_count_entries)
            settings.molecule_lengths.insert(settings.molecule_lengths.end(),
                                             static_cast<std::size_t>(entry.second), entry.first);
    } else {
        if (settings.length <= 0)
            throw std::runtime_error("--length must be positive");
        if (settings.chains <= 0)
            throw std::runtime_error("--chains must be positive");
        settings.total_repeats = 1LL * settings.length * settings.chains;
        settings.molecule_lengths.assign(static_cast<std::size_t>(settings.chains),
                                         settings.length);
    }
    if (settings.mps_weight_percent >= 0.0 && settings.mps_percent_explicit)
        throw std::runtime_error("--mps-wt and --mps-percent are mutually exclusive");

    double requested_mps_count = 0.0;
    if (settings.mps_weight_percent >= 0.0) {
        if (settings.mps_weight_percent > 100.0)
            throw std::runtime_error("--mps-wt must be between 0 and 100");
        const double fraction = settings.mps_weight_percent / 100.0;
        const double denominator =
            kMpsRepeatMass * (1.0 - fraction) + fraction * kDmsMass;
        requested_mps_count =
            denominator == 0.0 ? settings.total_repeats
                               : fraction * settings.total_repeats * kDmsMass / denominator;
    } else {
        if (settings.mps_monomer_percent < 0.0 ||
            settings.mps_monomer_percent > 100.0)
            throw std::runtime_error("--mps-percent must be between 0 and 100");
        requested_mps_count =
            settings.total_repeats * settings.mps_monomer_percent / 100.0;
    }
    settings.total_mps_repeats = clamp_value(
        std::llround(requested_mps_count), 0LL, settings.total_repeats);
    settings.base_mps_per_chain = static_cast<int>(
        settings.total_mps_repeats / settings.chains);
    settings.chains_with_extra_mps = static_cast<int>(
        settings.total_mps_repeats % settings.chains);
    settings.molecule_mps_counts.assign(
        static_cast<std::size_t>(settings.chains), settings.base_mps_per_chain);
    if (!settings.length_list_mode) {
        for (int i = 0; i < settings.chains_with_extra_mps; ++i)
            ++settings.molecule_mps_counts[static_cast<std::size_t>(i)];
    } else {
        // Largest remainders distribute a global composition across unequal chains.
        std::vector<std::pair<double, int>> remainders;
        remainders.reserve(static_cast<std::size_t>(settings.chains));
        long long allocated = 0;
        for (int i = 0; i < settings.chains; ++i) {
            const double quota = static_cast<double>(settings.molecule_lengths[static_cast<std::size_t>(i)]) *
                                 settings.total_mps_repeats / settings.total_repeats;
            const int base = static_cast<int>(std::floor(quota));
            settings.molecule_mps_counts[static_cast<std::size_t>(i)] = base;
            allocated += base;
            remainders.emplace_back(quota - base, i);
        }
        std::sort(remainders.begin(), remainders.end(),
                  [](const auto& a, const auto& b) {
                      if (a.first != b.first) return a.first > b.first;
                      return a.second < b.second;
                  });
        for (long long i = 0; i < settings.total_mps_repeats - allocated; ++i)
            ++settings.molecule_mps_counts[static_cast<std::size_t>(
                remainders[static_cast<std::size_t>(i)].second)];
    }
}

void validate(const Settings& settings) {
    if (settings.chains <= 0)
        throw std::runtime_error("--chains must be positive");
    if (settings.sequence != "random" &&
        settings.sequence != "alternating" &&
        settings.sequence != "block")
        throw std::runtime_error("--sequence must be random, alternating, or block");
    if (settings.density <= 0.0)
        throw std::runtime_error("--density must be positive");
    if (settings.target_density <= 0.0)
        throw std::runtime_error("--target-density must be positive");
    if (settings.target_density < settings.density)
        throw std::runtime_error("--target-density must be at least --density");
    if (settings.film_padding != -1.0 && settings.film_padding <= 0.0)
        throw std::runtime_error("--film-padding must be positive");
    if (settings.minimum_separation <= 0.0)
        throw std::runtime_error("--min-separation must be positive");
    if (settings.minimum_separation >= 15.0)
        throw std::runtime_error("--min-separation must be less than the 15 A pair cutoff");
    if (settings.output_explicit && settings.output.empty())
        throw std::runtime_error("--output cannot be empty");

    const long long maximum_atoms =
        settings.total_repeats + settings.total_mps_repeats;
    if (maximum_atoms > std::numeric_limits<int>::max())
        throw std::runtime_error("The requested system exceeds 32-bit LAMMPS atom IDs");
}

void derive_output_name(Settings& settings) {
    if (settings.output_explicit) return;
    std::ostringstream name;
    if (settings.total_mps_repeats == 0) {
        name << "data.Oil_PDMS";
    } else if (settings.total_mps_repeats == settings.total_repeats) {
        name << "data.Oil_PMPS";
    } else {
        name << "data.Oil_Copolymer";
    }
    if (settings.length_list_mode)
        name << "_R" << settings.total_repeats << "_M" << settings.chains;
    else
        name << "_N" << settings.length << "_M" << settings.chains;
    if (settings.total_mps_repeats != 0 &&
        settings.total_mps_repeats != settings.total_repeats) {
        if (settings.mps_weight_percent >= 0.0)
            name << "_MPS" << filename_number(settings.mps_weight_percent) << "wt";
        else
            name << "_MPS" << filename_number(settings.mps_monomer_percent) << "mon";
        name << '_' << settings.sequence;
    }
    settings.output = name.str();
}

double total_mass(const Settings& settings) {
    return (settings.total_repeats - settings.total_mps_repeats) * kDmsMass +
           settings.total_mps_repeats * kMpsRepeatMass;
}

Box calculate_box(const Settings& settings) {
    const double volume = total_mass(settings) / (kAvogadroScale * settings.density);
    const double length = std::cbrt(volume);
    return {length, length, length};
}

std::vector<bool> make_sequence(
    const Settings& settings,
    int length,
    int count,
    std::mt19937& rng
) {
    std::vector<bool> is_mps(static_cast<std::size_t>(length), false);
    if (count == 0) return is_mps;
    if (count == length) {
        std::fill(is_mps.begin(), is_mps.end(), true);
        return is_mps;
    }

    if (settings.sequence == "random") {
        std::vector<int> sites(static_cast<std::size_t>(length));
        for (int i = 0; i < length; ++i) sites[static_cast<std::size_t>(i)] = i;
        std::shuffle(sites.begin(), sites.end(), rng);
        for (int i = 0; i < count; ++i)
            is_mps[static_cast<std::size_t>(sites[static_cast<std::size_t>(i)])] = true;
    } else if (settings.sequence == "alternating") {
        for (int i = 0; i < length; ++i) {
            const int before = i * count / length;
            const int after = (i + 1) * count / length;
            if (after > before) is_mps[static_cast<std::size_t>(i)] = true;
        }
    } else {
        const int first = (length - count) / 2;
        for (int i = first; i < first + count; ++i)
            is_mps[static_cast<std::size_t>(i)] = true;
    }
    return is_mps;
}

double backbone_bond_length(bool left_mps, bool right_mps) {
    return left_mps && right_mps
        ? kMpsBackboneBondLength
        : kDmsBondLength;
}

double backbone_angle_degrees(bool left_mps, bool center_mps, bool right_mps) {
    return left_mps && center_mps && right_mps
        ? kMpsBackboneAngleDegrees
        : kDmsAngleDegrees;
}

std::pair<Vec3, Vec3> perpendicular_basis(const Vec3& axis) {
    const Vec3 unit = normalized(axis);
    const Vec3 reference =
        std::fabs(unit.z) < 0.85 ? Vec3{0.0, 0.0, 1.0}
                                : Vec3{0.0, 1.0, 0.0};
    const Vec3 first = normalized(cross(reference, unit));
    const Vec3 second = normalized(cross(unit, first));
    return {first, second};
}

Vec3 random_unit_vector(std::mt19937& rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double z = 2.0 * unit(rng) - 1.0;
    const double phi = 2.0 * kPi * unit(rng);
    const double radius = std::sqrt(std::max(0.0, 1.0 - z * z));
    return {radius * std::cos(phi), radius * std::sin(phi), z};
}

int topological_distance(const LocalSite& a, const LocalSite& b) {
    const int backbone_steps = std::abs(a.monomer - b.monomer);
    return backbone_steps +
           static_cast<int>(a.pendant) +
           static_cast<int>(b.pendant);
}

bool local_pair_allowed(const LocalSite& a, const LocalSite& b) {
    const int path_length = topological_distance(a, b);
    if (path_length <= 2) return true;
    const double minimum_distance = path_length == 3 ? 4.0 : 4.5;
    return norm2(a.position - b.position) >=
           minimum_distance * minimum_distance;
}

bool candidate_sites_allowed(
    const std::vector<LocalSite>& existing,
    const std::vector<LocalSite>& candidates
) {
    for (const LocalSite& candidate : candidates) {
        for (const LocalSite& site : existing) {
            if (!local_pair_allowed(candidate, site)) return false;
        }
    }
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        for (std::size_t j = i + 1; j < candidates.size(); ++j) {
            if (!local_pair_allowed(candidates[i], candidates[j])) return false;
        }
    }
    return true;
}

Vec3 endpoint_pendant_position(
    const Vec3& center,
    const Vec3& neighbor,
    double phi
) {
    const Vec3 toward_neighbor = normalized(neighbor - center);
    const auto basis = perpendicular_basis(toward_neighbor);
    const double theta = radians(kMpsPendantAngleDegrees);
    const Vec3 direction =
        std::cos(theta) * toward_neighbor +
        std::sin(theta) *
            (std::cos(phi) * basis.first + std::sin(phi) * basis.second);
    return center + kMpsPendantBondLength * normalized(direction);
}

Vec3 interior_pendant_position(
    const Vec3& left,
    const Vec3& center,
    const Vec3& right,
    double sign
) {
    const Vec3 toward_left = normalized(left - center);
    const Vec3 toward_right = normalized(right - center);
    const Vec3 bisector = normalized(toward_left + toward_right);
    Vec3 normal = cross(toward_left, toward_right);
    if (norm(normal) < 1.0e-10)
        normal = perpendicular_basis(bisector).first;
    else
        normal = normalized(normal);
    const double theta = radians(kMpsPendantAngleDegrees);
    const double projection_denominator = dot(bisector, toward_left);
    double bisector_scale = std::cos(theta) / projection_denominator;
    bisector_scale = clamp_value(bisector_scale, -1.0, 1.0);
    const double normal_scale =
        std::sqrt(std::max(0.0, 1.0 - bisector_scale * bisector_scale));
    const Vec3 direction =
        bisector_scale * bisector + sign * normal_scale * normal;
    return center + kMpsPendantBondLength * normalized(direction);
}

LocalMolecule build_local_molecule(
    const std::vector<bool>& is_mps,
    std::mt19937& rng
) {
    const int length = static_cast<int>(is_mps.size());
    LocalMolecule molecule;
    molecule.backbone.resize(static_cast<std::size_t>(length));
    molecule.pendants.resize(static_cast<std::size_t>(length));
    std::vector<LocalSite> placed;
    placed.reserve(static_cast<std::size_t>(
        length + std::count(is_mps.begin(), is_mps.end(), true)));
    std::uniform_real_distribution<double> azimuth(0.0, 2.0 * kPi);
    const int side_parity = std::uniform_int_distribution<int>(0, 1)(rng);

    placed.push_back({molecule.backbone[0], 0, false});
    if (length == 1) {
        if (is_mps[0]) {
            molecule.pendants[0] =
                kMpsPendantBondLength * random_unit_vector(rng);
        }
    } else {
        const double first_length = backbone_bond_length(is_mps[0], is_mps[1]);
        molecule.backbone[1] = {first_length, 0.0, 0.0};
        placed.push_back({molecule.backbone[1], 1, false});
        if (is_mps[0]) {
            molecule.pendants[0] = endpoint_pendant_position(
                molecule.backbone[0],
                molecule.backbone[1],
                azimuth(rng));
            placed.push_back({molecule.pendants[0], 0, true});
        }

        for (int i = 2; i < length; ++i) {
            bool accepted = false;
            for (int attempt = 0; attempt < 200 && !accepted; ++attempt) {
                const Vec3 previous_bond = normalized(
                    molecule.backbone[static_cast<std::size_t>(i - 1)] -
                    molecule.backbone[static_cast<std::size_t>(i - 2)]);
                const auto basis = perpendicular_basis(previous_bond);
                const double theta = radians(backbone_angle_degrees(
                    is_mps[static_cast<std::size_t>(i - 2)],
                    is_mps[static_cast<std::size_t>(i - 1)],
                    is_mps[static_cast<std::size_t>(i)]));
                const double alpha = kPi - theta;
                const double phi = azimuth(rng);
                const Vec3 direction =
                    std::cos(alpha) * previous_bond +
                    std::sin(alpha) *
                        (std::cos(phi) * basis.first +
                         std::sin(phi) * basis.second);
                const double bond_length = backbone_bond_length(
                    is_mps[static_cast<std::size_t>(i - 1)],
                    is_mps[static_cast<std::size_t>(i)]);
                const Vec3 candidate_backbone =
                    molecule.backbone[static_cast<std::size_t>(i - 1)] +
                    bond_length * normalized(direction);
                std::vector<LocalSite> candidates;
                candidates.push_back({candidate_backbone, i, false});
                Vec3 candidate_pendant;
                if (is_mps[static_cast<std::size_t>(i - 1)]) {
                    const double sign =
                        ((i - 1 + side_parity) % 2 == 0) ? 1.0 : -1.0;
                    candidate_pendant = interior_pendant_position(
                        molecule.backbone[static_cast<std::size_t>(i - 2)],
                        molecule.backbone[static_cast<std::size_t>(i - 1)],
                        candidate_backbone,
                        sign);
                    candidates.push_back({candidate_pendant, i - 1, true});
                }
                if (!candidate_sites_allowed(placed, candidates)) continue;
                accepted = true;
                molecule.backbone[static_cast<std::size_t>(i)] =
                    candidate_backbone;
                placed.insert(placed.end(), candidates.begin(), candidates.end());
                if (is_mps[static_cast<std::size_t>(i - 1)])
                    molecule.pendants[static_cast<std::size_t>(i - 1)] =
                        candidate_pendant;
            }
            if (!accepted) {
                molecule.backbone.clear();
                molecule.pendants.clear();
                return molecule;
            }
        }

        if (is_mps[static_cast<std::size_t>(length - 1)]) {
            bool accepted = false;
            for (int attempt = 0; attempt < 200 && !accepted; ++attempt) {
                const Vec3 candidate = endpoint_pendant_position(
                    molecule.backbone[static_cast<std::size_t>(length - 1)],
                    molecule.backbone[static_cast<std::size_t>(length - 2)],
                    azimuth(rng));
                const std::vector<LocalSite> candidates{
                    {candidate, length - 1, true}
                };
                if (!candidate_sites_allowed(placed, candidates)) continue;
                molecule.pendants[static_cast<std::size_t>(length - 1)] =
                    candidate;
                placed.push_back(candidates.front());
                accepted = true;
            }
            if (!accepted) {
                molecule.backbone.clear();
                molecule.pendants.clear();
                return molecule;
            }
        }
    }

    Vec3 center;
    int atom_count = 0;
    for (int i = 0; i < length; ++i) {
        center = center + molecule.backbone[static_cast<std::size_t>(i)];
        ++atom_count;
        if (is_mps[static_cast<std::size_t>(i)]) {
            center = center + molecule.pendants[static_cast<std::size_t>(i)];
            ++atom_count;
        }
    }
    center = center / static_cast<double>(atom_count);
    for (int i = 0; i < length; ++i) {
        molecule.backbone[static_cast<std::size_t>(i)] =
            molecule.backbone[static_cast<std::size_t>(i)] - center;
        if (is_mps[static_cast<std::size_t>(i)])
            molecule.pendants[static_cast<std::size_t>(i)] =
                molecule.pendants[static_cast<std::size_t>(i)] - center;
    }
    return molecule;
}

std::vector<LocalSite> local_sites(
    const LocalMolecule& molecule,
    const std::vector<bool>& is_mps
) {
    std::vector<LocalSite> sites;
    sites.reserve(molecule.backbone.size() +
                  static_cast<std::size_t>(std::count(is_mps.begin(), is_mps.end(), true)));
    for (std::size_t i = 0; i < molecule.backbone.size(); ++i)
        sites.push_back({molecule.backbone[i], static_cast<int>(i), false});
    for (std::size_t i = 0; i < molecule.pendants.size(); ++i) {
        if (is_mps[i])
            sites.push_back({molecule.pendants[i], static_cast<int>(i), true});
    }
    return sites;
}

bool internally_valid(const std::vector<LocalSite>& sites) {
    for (std::size_t i = 0; i < sites.size(); ++i) {
        for (std::size_t j = i + 1; j < sites.size(); ++j) {
            if (!local_pair_allowed(sites[i], sites[j])) return false;
        }
    }
    return true;
}

std::array<double, 9> random_rotation(std::mt19937& rng) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double u1 = unit(rng);
    const double u2 = unit(rng);
    const double u3 = unit(rng);
    const double qx = std::sqrt(1.0 - u1) * std::sin(2.0 * kPi * u2);
    const double qy = std::sqrt(1.0 - u1) * std::cos(2.0 * kPi * u2);
    const double qz = std::sqrt(u1) * std::sin(2.0 * kPi * u3);
    const double qw = std::sqrt(u1) * std::cos(2.0 * kPi * u3);
    return {
        1.0 - 2.0 * (qy * qy + qz * qz),
        2.0 * (qx * qy - qz * qw),
        2.0 * (qx * qz + qy * qw),
        2.0 * (qx * qy + qz * qw),
        1.0 - 2.0 * (qx * qx + qz * qz),
        2.0 * (qy * qz - qx * qw),
        2.0 * (qx * qz - qy * qw),
        2.0 * (qy * qz + qx * qw),
        1.0 - 2.0 * (qx * qx + qy * qy)
    };
}

Vec3 rotate(const std::array<double, 9>& matrix, const Vec3& v) {
    return {
        matrix[0] * v.x + matrix[1] * v.y + matrix[2] * v.z,
        matrix[3] * v.x + matrix[4] * v.y + matrix[5] * v.z,
        matrix[6] * v.x + matrix[7] * v.y + matrix[8] * v.z
    };
}

Vec3 minimum_image(Vec3 delta, const Box& box, bool film) {
    delta.x -= std::round(delta.x / box.lx) * box.lx;
    delta.y -= std::round(delta.y / box.ly) * box.ly;
    if (!film) delta.z -= std::round(delta.z / box.lz) * box.lz;
    return delta;
}

class PeriodicCellList {
public:
    PeriodicCellList(const Box& box, double cutoff, bool film)
        : box_(box), film_(film),
          cutoff2_(cutoff * cutoff),
          cell_counts_{{std::max(1, static_cast<int>(std::floor(box.lx / cutoff))),
                        std::max(1, static_cast<int>(std::floor(box.ly / cutoff))),
                        std::max(1, static_cast<int>(std::floor(box.lz / cutoff))) }},
          cell_widths_{{box.lx / cell_counts_[0], box.ly / cell_counts_[1],
                        box.lz / cell_counts_[2]}} {}

    bool overlaps(const std::vector<Vec3>& positions) const {
        for (const Vec3& position : positions) {
            const auto cell = cell_coordinates(position);
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        const int z = cell[2] + dz;
                        if (film_ && (z < 0 || z >= cell_counts_[2])) continue;
                        const long long key = cell_key(
                            wrapped_index(cell[0] + dx, 0),
                            wrapped_index(cell[1] + dy, 1),
                            film_ ? z : wrapped_index(z, 2));
                        const auto found = cells_.find(key);
                        if (found == cells_.end()) continue;
                        for (const Vec3& existing : found->second) {
                            if (norm2(minimum_image(position - existing, box_, film_)) <
                                cutoff2_)
                                return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    void insert(const std::vector<Vec3>& positions) {
        for (const Vec3& position : positions) {
            const auto cell = cell_coordinates(position);
            cells_[cell_key(cell[0], cell[1], cell[2])].push_back(position);
        }
    }

private:
    std::array<int, 3> cell_coordinates(const Vec3& position) const {
        auto index = [&](double value, double length, int axis) {
            int result = static_cast<int>(std::floor((value + 0.5 * length) /
                                                     cell_widths_[axis]));
            return clamp_value(result, 0, cell_counts_[axis] - 1);
        };
        return {index(position.x, box_.lx, 0), index(position.y, box_.ly, 1),
                index(position.z, box_.lz, 2)};
    }

    int wrapped_index(int index, int axis) const {
        index %= cell_counts_[axis];
        if (index < 0) index += cell_counts_[axis];
        return index;
    }

    long long cell_key(int x, int y, int z) const {
        return (static_cast<long long>(x) * cell_counts_[1] + y) *
               cell_counts_[2] + z;
    }

    Box box_;
    bool film_;
    double cutoff2_;
    std::array<int, 3> cell_counts_;
    std::array<double, 3> cell_widths_;
    std::unordered_map<long long, std::vector<Vec3>> cells_;
};

int backbone_atom_type(bool is_mps) {
    return is_mps ? 4 : 1;
}

void add_topology(
    System& system,
    const std::vector<bool>& is_mps,
    const std::vector<int>& backbone_ids,
    const std::vector<int>& pendant_ids
) {
    const int length = static_cast<int>(is_mps.size());
    for (int i = 0; i + 1 < length; ++i) {
        const int type =
            is_mps[static_cast<std::size_t>(i)] &&
            is_mps[static_cast<std::size_t>(i + 1)] ? 2 : 1;
        system.bonds.push_back({
            static_cast<int>(system.bonds.size()) + 1,
            type,
            backbone_ids[static_cast<std::size_t>(i)],
            backbone_ids[static_cast<std::size_t>(i + 1)]
        });
    }
    for (int i = 0; i < length; ++i) {
        if (!is_mps[static_cast<std::size_t>(i)]) continue;
        system.bonds.push_back({
            static_cast<int>(system.bonds.size()) + 1,
            3,
            backbone_ids[static_cast<std::size_t>(i)],
            pendant_ids[static_cast<std::size_t>(i)]
        });
    }

    for (int i = 1; i + 1 < length; ++i) {
        const int type =
            is_mps[static_cast<std::size_t>(i - 1)] &&
            is_mps[static_cast<std::size_t>(i)] &&
            is_mps[static_cast<std::size_t>(i + 1)] ? 2 : 1;
        system.angles.push_back({
            static_cast<int>(system.angles.size()) + 1,
            type,
            backbone_ids[static_cast<std::size_t>(i - 1)],
            backbone_ids[static_cast<std::size_t>(i)],
            backbone_ids[static_cast<std::size_t>(i + 1)]
        });
    }
    for (int i = 0; i < length; ++i) {
        if (!is_mps[static_cast<std::size_t>(i)]) continue;
        if (i > 0) {
            system.angles.push_back({
                static_cast<int>(system.angles.size()) + 1,
                3,
                backbone_ids[static_cast<std::size_t>(i - 1)],
                backbone_ids[static_cast<std::size_t>(i)],
                pendant_ids[static_cast<std::size_t>(i)]
            });
        }
        if (i + 1 < length) {
            system.angles.push_back({
                static_cast<int>(system.angles.size()) + 1,
                3,
                pendant_ids[static_cast<std::size_t>(i)],
                backbone_ids[static_cast<std::size_t>(i)],
                backbone_ids[static_cast<std::size_t>(i + 1)]
            });
        }
    }

    for (int i = 0; i + 1 < length; ++i) {
        struct EndAtom {
            int id;
            bool pendant;
            bool mps_backbone;
        };
        std::vector<EndAtom> left;
        std::vector<EndAtom> right;
        if (i > 0) {
            left.push_back({
                backbone_ids[static_cast<std::size_t>(i - 1)],
                false,
                is_mps[static_cast<std::size_t>(i - 1)]
            });
        }
        if (is_mps[static_cast<std::size_t>(i)]) {
            left.push_back({
                pendant_ids[static_cast<std::size_t>(i)],
                true,
                false
            });
        }
        if (i + 2 < length) {
            right.push_back({
                backbone_ids[static_cast<std::size_t>(i + 2)],
                false,
                is_mps[static_cast<std::size_t>(i + 2)]
            });
        }
        if (is_mps[static_cast<std::size_t>(i + 1)]) {
            right.push_back({
                pendant_ids[static_cast<std::size_t>(i + 1)],
                true,
                false
            });
        }

        for (const EndAtom& a : left) {
            for (const EndAtom& d : right) {
                const int pendant_count =
                    static_cast<int>(a.pendant) + static_cast<int>(d.pendant);
                int type = 1;
                if (pendant_count == 2) {
                    type = 4;
                } else if (pendant_count == 1) {
                    type = 3;
                } else {
                    const bool all_mps =
                        a.mps_backbone &&
                        is_mps[static_cast<std::size_t>(i)] &&
                        is_mps[static_cast<std::size_t>(i + 1)] &&
                        d.mps_backbone;
                    type = all_mps ? 2 : 1;
                }
                system.dihedrals.push_back({
                    static_cast<int>(system.dihedrals.size()) + 1,
                    type,
                    a.id,
                    backbone_ids[static_cast<std::size_t>(i)],
                    backbone_ids[static_cast<std::size_t>(i + 1)],
                    d.id
                });
            }
        }
    }
}

System generate_system(
    const Settings& settings,
    const Box& box
) {
    System system;
    const std::size_t expected_atoms = static_cast<std::size_t>(
        settings.total_repeats + settings.total_mps_repeats);
    system.atoms.reserve(expected_atoms);

    std::mt19937 rng(settings.seed);
    std::vector<std::pair<int, int>> molecules;
    molecules.reserve(static_cast<std::size_t>(settings.chains));
    for (int i = 0; i < settings.chains; ++i)
        molecules.emplace_back(settings.molecule_lengths[static_cast<std::size_t>(i)],
                               settings.molecule_mps_counts[static_cast<std::size_t>(i)]);
    if (settings.length_list_mode || settings.chains_with_extra_mps > 0)
        std::shuffle(molecules.begin(), molecules.end(), rng);
    PeriodicCellList cell_list(box, settings.minimum_separation, false);
    const int nx = static_cast<int>(std::ceil(std::cbrt(
        static_cast<double>(settings.chains) * box.lx / box.lz)));
    const int ny = nx;
    const int nz = (settings.chains + nx * ny - 1) / (nx * ny);
    const Vec3 cell_spacing = {
        box.lx / nx,
        box.ly / ny,
        box.lz / nz
    };
    std::uniform_real_distribution<double> jitter(-0.15, 0.15);

    for (int molecule_index = 0; molecule_index < settings.chains; ++molecule_index) {
        const int length = molecules[static_cast<std::size_t>(molecule_index)].first;
        const std::vector<bool> is_mps = make_sequence(
            settings, length,
            molecules[static_cast<std::size_t>(molecule_index)].second, rng);
        const int gx = molecule_index % nx;
        const int gy = (molecule_index / nx) % ny;
        const int gz = molecule_index / (nx * ny);
        const Vec3 base_anchor = {
            -0.5 * box.lx + (gx + 0.5) * cell_spacing.x,
            -0.5 * box.ly + (gy + 0.5) * cell_spacing.y,
            -0.5 * box.lz + (gz + 0.5) * cell_spacing.z
        };

        std::vector<Vec3> accepted_positions;
        bool accepted = false;
        for (int attempt = 0; attempt < 400 && !accepted; ++attempt) {
            LocalMolecule local = build_local_molecule(is_mps, rng);
            if (local.backbone.empty()) continue;
            const std::vector<LocalSite> sites = local_sites(local, is_mps);
            if (!internally_valid(sites)) continue;
            const auto rotation = random_rotation(rng);

            std::vector<Vec3> rotated;
            rotated.reserve(sites.size());
            Vec3 lower = {
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()
            };
            Vec3 upper = {
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest()
            };
            for (const LocalSite& site : sites) {
                const Vec3 position = rotate(rotation, site.position);
                rotated.push_back(position);
                lower.x = std::min(lower.x, position.x);
                lower.y = std::min(lower.y, position.y);
                lower.z = std::min(lower.z, position.z);
                upper.x = std::max(upper.x, position.x);
                upper.y = std::max(upper.y, position.y);
                upper.z = std::max(upper.z, position.z);
            }

            const Vec3 anchor_lower = {
                -0.5 * box.lx + kBoundaryClearance - lower.x,
                -0.5 * box.ly + kBoundaryClearance - lower.y,
                -0.5 * box.lz + kBoundaryClearance - lower.z
            };
            const Vec3 anchor_upper = {
                0.5 * box.lx - kBoundaryClearance - upper.x,
                0.5 * box.ly - kBoundaryClearance - upper.y,
                0.5 * box.lz - kBoundaryClearance - upper.z
            };
            if (anchor_lower.x > anchor_upper.x ||
                anchor_lower.y > anchor_upper.y ||
                anchor_lower.z > anchor_upper.z)
                continue;

            const Vec3 requested_anchor = {
                base_anchor.x + jitter(rng) * cell_spacing.x,
                base_anchor.y + jitter(rng) * cell_spacing.y,
                base_anchor.z + jitter(rng) * cell_spacing.z
            };
            const Vec3 anchor = {
                clamp_value(requested_anchor.x, anchor_lower.x, anchor_upper.x),
                clamp_value(requested_anchor.y, anchor_lower.y, anchor_upper.y),
                clamp_value(requested_anchor.z, anchor_lower.z, anchor_upper.z)
            };
            std::vector<Vec3> positions;
            positions.reserve(rotated.size());
            for (const Vec3& position : rotated)
                positions.push_back(anchor + position);

            if (cell_list.overlaps(positions)) continue;
            accepted = true;
            accepted_positions = std::move(positions);
        }
        if (!accepted) {
            std::ostringstream message;
            message << "Could not place molecule " << molecule_index + 1
                    << " wholly inside the bulk box without overlap. Lower --density, "
                       "--min-separation, or --length.";
            throw std::runtime_error(message.str());
        }

        cell_list.insert(accepted_positions);
        std::vector<int> backbone_ids(static_cast<std::size_t>(length));
        std::vector<int> pendant_ids(static_cast<std::size_t>(length), 0);
        std::size_t site_index = 0;
        for (int i = 0; i < length; ++i, ++site_index) {
            const Vec3& position = accepted_positions[site_index];
            const int id = static_cast<int>(system.atoms.size()) + 1;
            backbone_ids[static_cast<std::size_t>(i)] = id;
            system.atoms.push_back({
                id,
                molecule_index + 1,
                backbone_atom_type(is_mps[static_cast<std::size_t>(i)]),
                0.0,
                position,
                0,
                0,
                0
            });
        }
        for (int i = 0; i < length; ++i) {
            if (!is_mps[static_cast<std::size_t>(i)]) continue;
            const Vec3& position = accepted_positions[site_index++];
            const int id = static_cast<int>(system.atoms.size()) + 1;
            pendant_ids[static_cast<std::size_t>(i)] = id;
            system.atoms.push_back({
                id,
                molecule_index + 1,
                5,
                0.0,
                position,
                0,
                0,
                0
            });
        }
        add_topology(system, is_mps, backbone_ids, pendant_ids);
    }

    if (system.atoms.size() != expected_atoms)
        throw std::logic_error("Generated atom count does not match the requested composition");
    return system;
}

OutputFiles output_files(const Settings& settings) {
    OutputFiles files;
    const std::filesystem::path requested(settings.output);
    files.data_basename = requested.filename().string();
    if (files.data_basename.empty())
        throw std::runtime_error("Output filename cannot end with a separator");
    files.case_name =
        files.data_basename.rfind("data.", 0) == 0
            ? files.data_basename.substr(5)
            : files.data_basename;
    if (files.case_name.empty() || files.case_name == "." || files.case_name == "..")
        throw std::runtime_error("Cannot derive a case name from --output");
    const std::filesystem::path directory = requested.has_parent_path()
        ? requested.parent_path() : std::filesystem::path(".");
    files.directory = directory.string();
    files.input_basename = "in." + files.case_name;
    files.submit_basename = "submit." + files.case_name + ".sh";
    files.film_input_basename = "in." + files.case_name + ".film";
    files.film_submit_basename = "submit." + files.case_name + ".film.sh";
    files.pair_submit_basename = "submit." + files.case_name + ".pair.sh";
    files.info_basename = files.case_name + ".info";
    files.film_equilibration_energy_basename =
        "energy." + files.case_name + ".film_eq.dat";
    files.film_energy_basename = "energy." + files.case_name + ".film.dat";
    files.bulk_equilibrated_data_basename = "data." + files.case_name + ".npt_eq";
    files.film_initial_data_basename = "data." + files.case_name + ".film_initial";
    files.film_equilibrated_data_basename = "data." + files.case_name + ".film_eq";
    files.film_final_data_basename = "data." + files.case_name + ".film_final";
    files.unwrapped_dump_basename = "unwrap." + files.case_name + ".lammpstrj";
    auto in_directory = [&](const std::string& basename) {
        return (directory / basename).string();
    };
    files.data = in_directory(files.data_basename);
    files.input = in_directory(files.input_basename);
    files.submit = in_directory(files.submit_basename);
    files.film_input = in_directory(files.film_input_basename);
    files.film_submit = in_directory(files.film_submit_basename);
    files.pair_submit = in_directory(files.pair_submit_basename);
    files.info = in_directory(files.info_basename);
    return files;
}

void create_output_directory(const OutputFiles& files) {
    std::error_code error;
    std::filesystem::create_directories(files.directory, error);
    if (error)
        throw std::runtime_error("Cannot create output directory " + files.directory +
                                 ": " + error.message());
}

void write_data(
    const OutputFiles& files,
    const System& system,
    const Box& box
) {
    std::ofstream out(files.data);
    if (!out) throw std::runtime_error("Cannot open output data file: " + files.data);
    out << "LAMMPS data file for standalone PDMS/PMPS silicone oil\n\n"
        << system.atoms.size() << " atoms\n"
        << "5 atom types\n"
        << system.bonds.size() << " bonds\n"
        << "3 bond types\n"
        << system.angles.size() << " angles\n"
        << "3 angle types\n"
        << system.dihedrals.size() << " dihedrals\n"
        << "4 dihedral types\n\n"
        << std::fixed << std::setprecision(8)
        << -0.5 * box.lx << ' ' << 0.5 * box.lx << " xlo xhi\n"
        << -0.5 * box.ly << ' ' << 0.5 * box.ly << " ylo yhi\n"
        << -0.5 * box.lz << ' ' << 0.5 * box.lz << " zlo zhi\n\n"
        << "Masses\n\n"
        << "1 " << kDmsMass << "\n"
        << "2 " << kDmsMass << "\n"
        << "3 " << kDmsMass << "\n"
        << "4 " << kMpsBackboneMass << "\n"
        << "5 " << kMpsPendantMass << "\n\n"
        << "Atoms # full\n\n";
    for (const Atom& atom : system.atoms) {
        out << atom.id << ' ' << atom.molecule << ' ' << atom.type << ' '
            << atom.charge << ' '
            << atom.position.x << ' ' << atom.position.y << ' ' << atom.position.z << ' '
            << atom.ix << ' ' << atom.iy << ' ' << atom.iz << '\n';
    }
    out << "\nBonds\n\n";
    for (const Bond& bond : system.bonds)
        out << bond.id << ' ' << bond.type << ' ' << bond.a << ' ' << bond.b << '\n';
    out << "\nAngles\n\n";
    for (const Angle& angle : system.angles)
        out << angle.id << ' ' << angle.type << ' '
            << angle.a << ' ' << angle.b << ' ' << angle.c << '\n';
    out << "\nDihedrals\n\n";
    for (const Dihedral& dihedral : system.dihedrals)
        out << dihedral.id << ' ' << dihedral.type << ' '
            << dihedral.a << ' ' << dihedral.b << ' '
            << dihedral.c << ' ' << dihedral.d << '\n';
    if (!out) throw std::runtime_error("Failed while writing data file: " + files.data);
}

PairParameters dms_pair_parameters(double temperature) {
    if (std::fabs(temperature - 300.0) < 1.0e-12)
        return {kDmsEpsilon, kDmsSigma};
    const double shifted_temperature = temperature - 186.04682;
    const double exponential_argument = 0.00758 * shifted_temperature;
    const double denominator = 1.0 + std::exp(exponential_argument);
    const double epsilon =
        (4.77795 / denominator + 1.47169) * 0.350646;
    const double sigma =
        (7.86548e-05 * temperature + 1.27856) * 4.95013;
    return {epsilon, sigma};
}

PairParameters pair_parameters(int type_i, int type_j, double temperature) {
    if (type_i > type_j) std::swap(type_i, type_j);
    const PairParameters dms = dms_pair_parameters(temperature);
    const double epsilon_scale = dms.epsilon / kDmsEpsilon;
    const double sigma_scale = dms.sigma / kDmsSigma;
    const PairParameters mps_backbone = {
        kMpsBackboneEpsilon * epsilon_scale,
        kMpsBackboneSigma * sigma_scale
    };
    const PairParameters mps_pendant = {
        kMpsPendantEpsilon * epsilon_scale,
        kMpsPendantSigma * sigma_scale
    };
    const PairParameters mps_backbone_pendant = {
        kMpsBackbonePendantEpsilon * epsilon_scale,
        kMpsBackbonePendantSigma * sigma_scale
    };
    if (type_i <= 3 && type_j <= 3)
        return dms;
    if (type_i <= 3 && type_j == 4)
        return {
            kDmsMpsEpsilonFactor *
                std::sqrt(dms.epsilon * mps_backbone.epsilon),
            0.5 * (dms.sigma + mps_backbone.sigma)
        };
    if (type_i <= 3 && type_j == 5)
        return {
            kDmsMpsEpsilonFactor *
                std::sqrt(dms.epsilon * mps_pendant.epsilon),
            0.5 * (dms.sigma + mps_pendant.sigma)
        };
    if (type_i == 4 && type_j == 4)
        return mps_backbone;
    if (type_i == 4 && type_j == 5)
        return mps_backbone_pendant;
    if (type_i == 5 && type_j == 5)
        return mps_pendant;
    throw std::logic_error("Unsupported atom-type pair");
}

double repulsive_cutoff(const PairParameters& pair) {
    return pair.sigma * std::pow(2.0, 1.0 / 6.0);
}

double maximum_repulsive_cutoff(double temperature) {
    double maximum = 0.0;
    for (int i = 1; i <= 5; ++i) {
        for (int j = i; j <= 5; ++j) {
            maximum = std::max(
                maximum,
                repulsive_cutoff(pair_parameters(i, j, temperature)));
        }
    }
    return maximum;
}

void write_pair_matrix(
    std::ostream& out,
    double temperature,
    bool include_repulsive_cutoff
) {
    out << std::fixed << std::setprecision(9);
    for (int i = 1; i <= 5; ++i) {
        for (int j = i; j <= 5; ++j) {
            const PairParameters pair = pair_parameters(i, j, temperature);
            out << "pair_coeff      " << i << ' ' << j << ' '
                << pair.epsilon << ' ' << pair.sigma;
            if (include_repulsive_cutoff)
                out << ' ' << repulsive_cutoff(pair);
            out << '\n';
        }
    }
}

void write_input(
    const Settings& settings,
    const OutputFiles& files
) {
    std::ofstream out(files.input);
    if (!out) throw std::runtime_error("Cannot open LAMMPS input file: " + files.input);
    const double hot_temperature = 800.0;
    const double cold_temperature = 300.0;
    const double compression_scale =
        std::cbrt(settings.density / settings.target_density);
    const double hot_global_cutoff =
        maximum_repulsive_cutoff(hot_temperature);

    out << "# Generated by the standalone silicone-oil generator\n"
        << "# Types 1-3: DMS namespace; type 4: MPS backbone; type 5: MPS pendant\n"
        << "# Bulk stage; the dependent film job reads the 300 K NPT snapshot.\n\n"
        << "units           real\n"
        << "boundary        p p p\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     hybrid harmonic quartic\n"
        << "dihedral_style  nharmonic\n"
        << "special_bonds   lj 0 0 0.5\n"
        << std::fixed << std::setprecision(9)
        << "pair_style      lj/cut " << hot_global_cutoff << "\n"
        << "comm_modify     cutoff 15\n"
        << "neighbor        2.5 bin\n"
        << "neigh_modify    delay 5 every 1\n\n"
        << "read_data       " << files.data_basename << "\n\n"
        << "# Bonds: DMS-like backbone, MPS-MPS backbone, MPS pendant\n"
        << "bond_coeff      1 115.4086 2.801\n"
        << "bond_coeff      2 108.3835 2.8039\n"
        << "bond_coeff      3 232.8302 3.12497\n\n"
        << "# Backbone-only mixed terms use DMS; pendant-containing terms use MPS\n"
        << "angle_coeff     1 harmonic 64.62431 111.623\n"
        << "angle_coeff     2 quartic 110.566 64.3974 -139.5241 80.974\n"
        << "angle_coeff     3 quartic 110.746 -20.8906 23.3707 180.3228\n\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n"
        << "dihedral_coeff  2 8 1.3730 0.2686 0.4017 -1.7250 -0.7052 4.1390 0.2327 -2.2635\n"
        << "dihedral_coeff  3 8 2.3494 -1.5840 -1.6463 3.2133 4.1479 -2.8795 -2.2665 1.1811\n"
        << "dihedral_coeff  4 8 2.23125 0.24735 2.4327 -2.8832 -4.7124 7.17825 2.33495 -3.74\n\n"
        << "# 800 K repulsive matrix. Each cutoff is 2^(1/6)*sigma for that pair.\n"
        << "# PMPS values use the V22/V35 DMS 800K/300K epsilon and sigma ratios.\n"
        << "# DMS-MPS mixing: 0.579966 * geometric epsilon; arithmetic sigma.\n";
    write_pair_matrix(out, hot_temperature, true);
    out << "\n"
        << "timestep        " << kTimestepFs << "\n"
        << "thermo          1000\n"
        << "thermo_style    custom step temp density lx ly lz pxx pyy pzz "
        << "etotal epair ebond eangle edihed\n"
        << "restart         100000 restart." << files.case_name
        << ".1 restart." << files.case_name << ".2\n\n"
        << "minimize        1e-4 1e-6 100 1000\n"
        << "velocity        all create 800.0 " << settings.velocity_seed
        << " mom yes rot yes dist gaussian\n\n"
        << "# 800 K initial relaxation\n"
        << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
        << "run             1000000\n"
        << "write_data      data." << files.case_name << ".rep_800 nocoeff\n"
        << "unfix           integrate\n"
        << "\n"
        << "# Isotropic compression from the initial density to the target density\n"
        << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
        << "fix             compress all deform 1 "
        << "x scale " << compression_scale << ' '
        << "y scale " << compression_scale << ' '
        << "z scale " << compression_scale << " units box\n"
        << "run             1000000\n"
        << "unfix           compress\n\n"
        << "# Relax at the compressed dimensions\n"
        << "run             1000000\n\n";

    out << "# Extended 800 K equilibration at fixed dimensions\n"
        << "run             2000000\n"
        << "unfix           integrate\n"
        << "write_data      data." << files.case_name << ".eq_800 nocoeff\n\n"
        << "# Switch to the explicit 300 K attractive matrix\n"
        << "pair_style      lj/gromacs 12 15\n";
    write_pair_matrix(out, cold_temperature, false);
    out << "\n"
        << "# Cool from 800 K to 300 K under isotropic NPT\n"
        << "fix             integrate all npt temp 800.0 300.0 50.0 "
        << "iso 1.0 1.0 500.0\n"
        << "run             1000000\n"
        << "write_data      data." << files.case_name << ".300 nocoeff\n"
        << "unfix           integrate\n\n"
        << "# Final 300 K equilibration under isotropic NPT\n"
        << "fix             integrate all npt temp 300.0 300.0 50.0 "
        << "iso 1.0 1.0 500.0\n"
        << "run             " << kBulkFinalNptSteps << "\n"
        << "write_data      " << files.bulk_equilibrated_data_basename << " nocoeff\n"
        << "unfix           integrate\n";
    if (!out) throw std::runtime_error("Failed while writing input file: " + files.input);
}

void write_film_input(const Settings& settings, const OutputFiles& files) {
    std::ofstream out(files.film_input);
    if (!out)
        throw std::runtime_error("Cannot open film input file: " + files.film_input);
    const PairParameters wall = dms_pair_parameters(300.0);
    const double padding = film_padding_angstrom(settings);

    out << std::fixed << std::setprecision(9)
        << "# Derived film: read the equilibrated bulk, preserve its topology and x/y box.\n"
        << "# Material slab thickness starts near bulk Lz; cell Lz grows by 2*padding.\n"
        << "units           real\n"
        << "boundary        p p p\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     hybrid harmonic quartic\n"
        << "dihedral_style  nharmonic\n"
        << "special_bonds   lj 0 0 0.5\n"
        << "pair_style      lj/gromacs 12 15\n"
        << "comm_modify     cutoff 15\n"
        << "neighbor        2.5 bin\n"
        << "neigh_modify    delay 5 every 1\n"
        << "read_data       " << files.bulk_equilibrated_data_basename << "\n\n"
        << "bond_coeff      1 115.4086 2.801\n"
        << "bond_coeff      2 108.3835 2.8039\n"
        << "bond_coeff      3 232.8302 3.12497\n"
        << "angle_coeff     1 harmonic 64.62431 111.623\n"
        << "angle_coeff     2 quartic 110.566 64.3974 -139.5241 80.974\n"
        << "angle_coeff     3 quartic 110.746 -20.8906 23.3707 180.3228\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n"
        << "dihedral_coeff  2 8 1.3730 0.2686 0.4017 -1.7250 -0.7052 4.1390 0.2327 -2.2635\n"
        << "dihedral_coeff  3 8 2.3494 -1.5840 -1.6463 3.2133 4.1479 -2.8795 -2.2665 1.1811\n"
        << "dihedral_coeff  4 8 2.23125 0.24735 2.4327 -2.8832 -4.7124 7.17825 2.33495 -3.74\n";
    write_pair_matrix(out, 300.0, false);
    out << "\n"
        << "# Keep bonded molecules intact while changing z from periodic to fixed.\n"
        << "# This follows the LAMMPS bulk-to-slab reset/image + unwrapped-dump method.\n"
        << "reset_atoms     image all\n"
        << "compute         zu_atom all property/atom zu\n"
        << "compute         zu_min all reduce min c_zu_atom\n"
        << "compute         zu_max all reduce max c_zu_atom\n"
        << "run             0\n"
        << "# Use at least 20 A on each face, enlarged if an unwrapped\n"
        << "# chain would otherwise touch or cross a safety wall.\n"
        << "# max(x,y) is not a scalar function in LAMMPS equal-style variables.\n"
        << "variable        lower_pad equal zlo-c_zu_min+" << repulsive_cutoff(wall) << "\n"
        << "variable        upper_pad equal c_zu_max-zhi+" << repulsive_cutoff(wall) << "\n"
        << "variable        edge_pad equal "
           "0.5*(v_lower_pad+v_upper_pad+abs(v_lower_pad-v_upper_pad))\n"
        << "variable        needed_pad equal 0.5*(" << padding
        << "+v_edge_pad+abs(" << padding << "-v_edge_pad))\n"
        << "variable        pad_used equal ${needed_pad}\n"
        << "print           \"Film z padding per face: ${pad_used} A\"\n"
        << "uncompute       zu_max\n"
        << "uncompute       zu_min\n"
        << "uncompute       zu_atom\n"
        << "write_dump      all custom " << files.unwrapped_dump_basename
        << " id xu yu zu\n"
        << "change_box      all z delta -${pad_used} ${pad_used}"
        << " boundary p p f units box\n"
        << "read_dump       " << files.unwrapped_dump_basename
        << " 0 x y z box no replace yes\n"
        << "write_data      " << files.film_initial_data_basename << " nocoeff\n\n"
        << "# Edge safety walls remain active through film production.\n"
        << "fix             zlo_wall all wall/lj126 zlo EDGE "
        << wall.epsilon << ' ' << wall.sigma << ' '
        << repulsive_cutoff(wall) << " units box\n"
        << "fix             zhi_wall all wall/lj126 zhi EDGE "
        << wall.epsilon << ' ' << wall.sigma << ' '
        << repulsive_cutoff(wall) << " units box\n"
        << "fix_modify      zlo_wall virial no\n"
        << "fix_modify      zhi_wall virial no\n"
        << "timestep        " << kTimestepFs << "\n"
        << "thermo          1000\n"
        << "thermo_style    custom step temp pe density lx ly lz pxx pyy pzz\n"
        << "# One brief trajectory for checking the initial guarded stage.\n"
        << "dump            filmtraj all custom " << kFilmWallSteps << " dump." << files.case_name
        << ".film.lammpstrj id mol type x y z ix iy iz\n"
        << "dump_modify     filmtraj sort id\n"
        << "velocity        all create 300.0 " << settings.velocity_seed
        << " mom yes rot yes dist gaussian\n"
        << "fix             integrate all npt temp 300.0 300.0 50.0 "
           "x 1.0 1.0 500.0 y 1.0 1.0 500.0 couple xy\n"
        << "run             " << kFilmWallSteps << "\n"
        << "unfix           integrate\n"
        << "undump          filmtraj\n\n"
        << "# Safety walls are distant guards; check that their forces stay zero.\n"
        << "reset_timestep  0 time 0.0\n"
        << "thermo          100000\n"
        << "thermo_style    custom time temp pe pxx pyy pzz lx ly lz\n"
        << "thermo_modify   format float %.12g\n"
        << "# Short names keep each quoted output command on one physical line.\n"
        << "variable        t equal time\n"
        << "variable        T equal temp\n"
        << "variable        e equal pe\n"
        << "variable        x equal pxx\n"
        << "variable        y equal pyy\n"
        << "variable        z equal pzz\n"
        << "variable        a equal lx\n"
        << "variable        b equal ly\n"
        << "variable        c equal lz\n"
        << "variable        l equal f_zlo_wall[1]\n"
        << "variable        h equal f_zhi_wall[1]\n"
        << "fix             integrate all nvt temp 300.0 300.0 50.0\n"
        << "fix             equil_output all print " << kEnergySampleEverySteps
        << " \"$t $T $e $x $y $z $a $b $c $l $h\" file "
        << files.film_equilibration_energy_basename
        << " screen no title \"# time_fs temp_K pe_kcal_per_mol pxx_atm pyy_atm "
           "pzz_atm lx_A ly_A lz_A wall_lo_force wall_hi_force\"\n"
        << "run             " << kFilmRelaxSteps << "\n"
        << "unfix           equil_output\n"
        << "unfix           integrate\n"
        << "write_data      " << files.film_equilibrated_data_basename << " nocoeff\n"
        << "reset_timestep  0 time 0.0\n"
        << "fix             integrate all nvt temp 300.0 300.0 50.0\n"
        << "fix             energy_output all print " << kEnergySampleEverySteps
        << " \"$t $T $e $x $y $z $a $b $c $l $h\" file "
        << files.film_energy_basename
        << " screen no title \"# time_fs temp_K pe_kcal_per_mol pxx_atm pyy_atm "
           "pzz_atm lx_A ly_A lz_A wall_lo_force wall_hi_force\"\n"
        << "run             " << kFilmProductionSteps << "\n"
        << "unfix           energy_output\n"
        << "unfix           integrate\n"
        << "write_data      " << files.film_final_data_basename << " nocoeff\n";
    if (!out)
        throw std::runtime_error("Failed while writing film input file: " + files.film_input);
}

std::string sanitize_job_name(std::string name) {
    for (char& character : name) {
        const bool valid =
            (character >= 'a' && character <= 'z') ||
            (character >= 'A' && character <= 'Z') ||
            (character >= '0' && character <= '9') ||
            character == '_' || character == '-';
        if (!valid) character = '_';
    }
    if (name.size() > 100) name.resize(100);
    return name;
}

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char character : value) {
        if (character == '\'') quoted += "'\\''";
        else quoted.push_back(character);
    }
    return quoted + "'";
}

void write_submit_file(const std::string& path, const std::string& job_name,
                       const std::string& input_basename,
                       const std::string& output_basename,
                       const std::string& required_data = "") {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open Slurm file: " + path);
    out << "#!/bin/bash\n"
        << "#SBATCH --job-name=" << sanitize_job_name(job_name) << "\n"
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
        << "INPUT=" << shell_quote(input_basename) << "\n"
        << "OUTPUT=" << shell_quote(output_basename) << "\n\n";
    if (!required_data.empty())
        out << "if [[ ! -s " << shell_quote(required_data) << " ]]; then\n"
            << "  echo 'Missing equilibrated bulk data; run the bulk job first.' >&2\n"
            << "  exit 1\n"
            << "fi\n\n";
    out << "srun lmp -in \"$INPUT\" > \"$OUTPUT\"\n";
    if (!out) throw std::runtime_error("Failed while writing Slurm file: " + path);
}

void write_submit(const OutputFiles& files) {
    write_submit_file(files.submit, files.case_name + "_bulk",
                      files.input_basename, "out." + files.case_name + ".bulk");
    write_submit_file(files.film_submit, files.case_name + "_film",
                      files.film_input_basename, "out." + files.case_name + ".film",
                      files.bulk_equilibrated_data_basename);
    std::ofstream out(files.pair_submit);
    if (!out)
        throw std::runtime_error("Cannot open pair submit script: " + files.pair_submit);
    out << "#!/bin/bash\n"
        << "set -euo pipefail\n"
        << "cd -- \"$(dirname -- \"${BASH_SOURCE[0]}\")\"\n"
        << "bulk_job_id=$(sbatch --parsable " << shell_quote(files.submit_basename)
        << ")\n"
        << "bulk_job_id=${bulk_job_id%%;*}\n"
        << "film_job_id=$(sbatch --parsable --dependency=afterok:\"$bulk_job_id\" "
        << shell_quote(files.film_submit_basename) << ")\n"
        << "printf 'Bulk job: %s\\nFilm job: %s\\n' \"$bulk_job_id\" \"$film_job_id\"\n";
    if (!out)
        throw std::runtime_error("Failed while writing pair submit script: " +
                                 files.pair_submit);
}

std::array<long long, 6> atom_type_counts(const System& system) {
    std::array<long long, 6> counts{};
    for (const Atom& atom : system.atoms)
        ++counts[static_cast<std::size_t>(atom.type)];
    return counts;
}

template <typename Interaction>
std::vector<long long> interaction_type_counts(
    const std::vector<Interaction>& interactions,
    int type_count
) {
    std::vector<long long> counts(static_cast<std::size_t>(type_count + 1), 0);
    for (const Interaction& interaction : interactions)
        ++counts[static_cast<std::size_t>(interaction.type)];
    return counts;
}

std::string json_escape(const std::string& text) {
    std::string result;
    for (char character : text) {
        if (character == '\\' || character == '"') result.push_back('\\');
        result.push_back(character);
    }
    return result;
}

void write_json_pair_matrix(
    std::ostream& out,
    double temperature,
    bool include_repulsive_cutoff,
    const std::string& indent
) {
    out << "{\n";
    bool first = true;
    for (int i = 1; i <= 5; ++i) {
        for (int j = i; j <= 5; ++j) {
            if (!first) out << ",\n";
            const PairParameters pair = pair_parameters(i, j, temperature);
            out << indent << "  \"" << i << '-' << j << "\": {\"epsilon\": "
                << pair.epsilon << ", \"sigma\": " << pair.sigma;
            if (include_repulsive_cutoff)
                out << ", \"cutoff\": " << repulsive_cutoff(pair);
            out << '}';
            first = false;
        }
    }
    out << "\n" << indent << '}';
}

void write_info(
    const Settings& settings,
    const OutputFiles& files,
    const System& system,
    const Box& box
) {
    std::ofstream out(files.info);
    if (!out) throw std::runtime_error("Cannot open info file: " + files.info);
    const auto atoms = atom_type_counts(system);
    const auto bonds = interaction_type_counts(system.bonds, 3);
    const auto angles = interaction_type_counts(system.angles, 3);
    const auto dihedrals = interaction_type_counts(system.dihedrals, 4);
    const double realized_monomer_percent =
        100.0 * settings.total_mps_repeats / settings.total_repeats;
    const double realized_weight_percent =
        100.0 * settings.total_mps_repeats * kMpsRepeatMass /
        total_mass(settings);
    const bool uniform_length = std::all_of(
        settings.molecule_lengths.begin(), settings.molecule_lengths.end(),
        [&](int length) { return length == settings.molecule_lengths.front(); });
    const bool uniform_composition = std::all_of(
        settings.molecule_mps_counts.begin(), settings.molecule_mps_counts.end(),
        [&](int count) { return count == settings.molecule_mps_counts.front(); });
    std::map<int, int> composition_distribution;
    for (int count : settings.molecule_mps_counts)
        ++composition_distribution[count];
    const double compression_scale =
        std::cbrt(settings.density / settings.target_density);
    const double padding = film_padding_angstrom(settings);

    out << std::fixed << std::setprecision(8)
        << "{\n"
        << "  \"model\": \"standalone PDMS/PMPS silicone oil\",\n"
        << "  \"case_name\": \"" << json_escape(files.case_name) << "\",\n"
        << "  \"files\": {\n"
        << "    \"data\": \"" << json_escape(files.data_basename) << "\",\n"
        << "    \"lammps_input\": \"" << json_escape(files.input_basename) << "\",\n"
        << "    \"slurm_submit\": \"" << json_escape(files.submit_basename) << "\",\n"
        << "    \"film_lammps_input\": \""
        << json_escape(files.film_input_basename) << "\",\n"
        << "    \"film_slurm_submit\": \""
        << json_escape(files.film_submit_basename) << "\",\n"
        << "    \"pair_submit\": \""
        << json_escape(files.pair_submit_basename) << "\",\n"
        << "    \"bulk_equilibrated_data\": \""
        << json_escape(files.bulk_equilibrated_data_basename) << "\",\n"
        << "    \"film_equilibration_output\": \""
        << json_escape(files.film_equilibration_energy_basename) << "\",\n"
        << "    \"film_energy_output\": \""
        << json_escape(files.film_energy_basename) << "\",\n"
        << "    \"model_info\": \"" << json_escape(files.info_basename) << "\"\n"
        << "  },\n"
        << "  \"generator_input\": {\n"
        << "    \"config_file\": ";
    if (settings.config_file.empty()) out << "null";
    else out << '"' << json_escape(settings.config_file) << '"';
    out << ",\n"
        << "    \"precedence\": \"defaults < config file < command line\"\n"
        << "  },\n"
        << "  \"composition\": {\n"
        << "    \"chain_length\": ";
    if (uniform_length) out << settings.molecule_lengths.front();
    else out << "null";
    out << ",\n"
        << "    \"chain_length_distribution\": [";
    for (std::size_t i = 0; i < settings.chain_count_entries.size(); ++i) {
        if (i != 0) out << ", ";
        out << "{\"length\": " << settings.chain_count_entries[i].first
            << ", \"chains\": " << settings.chain_count_entries[i].second << '}';
    }
    if (!settings.length_list_mode)
        out << "{\"length\": " << settings.length
            << ", \"chains\": " << settings.chains << '}';
    out << "],\n"
        << "    \"chain_count\": " << settings.chains << ",\n"
        << "    \"total_repeats\": " << settings.total_repeats << ",\n"
        << "    \"dms_repeats_total\": "
        << settings.total_repeats - settings.total_mps_repeats << ",\n"
        << "    \"mps_repeats_total\": " << settings.total_mps_repeats << ",\n"
        << "    \"dms_repeats_per_chain\": ";
    if (uniform_length && uniform_composition)
        out << settings.molecule_lengths.front() - settings.molecule_mps_counts.front();
    else out << "null";
    out << ",\n"
        << "    \"mps_repeats_per_chain\": ";
    if (uniform_length && uniform_composition) out << settings.molecule_mps_counts.front();
    else out << "null";
    out << ",\n"
        << "    \"chain_composition_distribution\": [";
    bool first_composition = true;
    for (const auto& entry : composition_distribution) {
        if (!first_composition) out << ", ";
        out << "{\"mps_repeats\": " << entry.first
            << ", \"chains\": " << entry.second << '}';
        first_composition = false;
    }
    out << "],\n"
        << "    \"sequence\": \"" << settings.sequence << "\",\n"
        << "    \"requested_mps_monomer_percent\": ";
    if (settings.mps_weight_percent < 0.0)
        out << settings.mps_monomer_percent;
    else
        out << "null";
    out << ",\n"
        << "    \"requested_mps_weight_percent\": ";
    if (settings.mps_weight_percent >= 0.0)
        out << settings.mps_weight_percent;
    else
        out << "null";
    out << ",\n"
        << "    \"realized_mps_monomer_percent\": "
        << realized_monomer_percent << ",\n"
        << "    \"realized_mps_weight_percent\": "
        << realized_weight_percent << ",\n"
        << "    \"chain_mass_g_per_mol\": ";
    if (uniform_length && uniform_composition)
        out << total_mass(settings) / settings.chains;
    else out << "null";
    out << ",\n"
        << "    \"mean_chain_mass_g_per_mol\": "
        << total_mass(settings) / settings.chains << "\n"
        << "  },\n"
        << "  \"atom_types\": {\n"
        << "    \"1\": {\"name\": \"neutral DMS\", \"mass\": "
        << kDmsMass << ", \"count\": " << atoms[1] << "},\n"
        << "    \"2\": {\"name\": \"reserved reactive DMS\", \"mass\": "
        << kDmsMass << ", \"count\": " << atoms[2] << "},\n"
        << "    \"3\": {\"name\": \"reserved reactive DMS\", \"mass\": "
        << kDmsMass << ", \"count\": " << atoms[3] << "},\n"
        << "    \"4\": {\"name\": \"MPS backbone\", \"mass\": "
        << kMpsBackboneMass << ", \"count\": " << atoms[4] << "},\n"
        << "    \"5\": {\"name\": \"MPS pendant\", \"mass\": "
        << kMpsPendantMass << ", \"count\": " << atoms[5] << "}\n"
        << "  },\n"
        << "  \"topology_counts\": {\n"
        << "    \"atoms\": " << system.atoms.size() << ",\n"
        << "    \"bonds\": {\"total\": " << system.bonds.size()
        << ", \"type_1\": " << bonds[1] << ", \"type_2\": " << bonds[2]
        << ", \"type_3\": " << bonds[3] << "},\n"
        << "    \"angles\": {\"total\": " << system.angles.size()
        << ", \"type_1\": " << angles[1] << ", \"type_2\": " << angles[2]
        << ", \"type_3\": " << angles[3] << "},\n"
        << "    \"dihedrals\": {\"total\": " << system.dihedrals.size()
        << ", \"type_1\": " << dihedrals[1] << ", \"type_2\": " << dihedrals[2]
        << ", \"type_3\": " << dihedrals[3] << ", \"type_4\": "
        << dihedrals[4] << "}\n"
        << "  },\n"
        << "  \"initial_box\": {\n"
        << "    \"geometry\": \"bulk\",\n"
        << "    \"boundary\": \"p p p\",\n"
        << "    \"lx_angstrom\": " << box.lx << ",\n"
        << "    \"ly_angstrom\": " << box.ly << ",\n"
        << "    \"lz_angstrom\": " << box.lz << ",\n"
        << "    \"length_angstrom\": " << box.lx << ",\n"
        << "    \"density_g_cm3\": " << settings.density << ",\n"
        << "    \"minimum_inter_molecular_separation_angstrom\": "
        << settings.minimum_separation << ",\n"
        << "    \"whole_chains_inside_primary_box\": true,\n"
        << "    \"atom_image_flags\": [0, 0, 0]\n"
        << "  },\n"
        << "  \"random_seeds\": {\n"
        << "    \"generator\": " << settings.seed << ",\n"
        << "    \"velocity\": " << settings.velocity_seed << "\n"
        << "  },\n"
        << "  \"mixing\": {\n"
        << "    \"dms_mps_rule\": \"0.579966 * geometric epsilon and arithmetic sigma\",\n"
        << "    \"dms_mps_epsilon_factor\": " << kDmsMpsEpsilonFactor << ",\n"
        << "    \"pmps_temperature_scaling\": "
        << "\"300 K PMPS epsilon and sigma scaled by the V22/V35 DMS 800K/300K ratios\",\n"
        << "    \"all_pair_coefficients_written_explicitly\": true\n"
        << "  },\n"
        << "  \"pair_coefficients_800K\": ";
    write_json_pair_matrix(out, 800.0, true, "  ");
    out << ",\n"
        << "  \"pair_coefficients_300K\": ";
    write_json_pair_matrix(out, 300.0, false, "  ");
    out << ",\n"
        << "  \"simulation_template\": {\n"
        << "    \"initial_density_g_cm3\": " << settings.density << ",\n"
        << "    \"target_compressed_density_g_cm3\": "
        << settings.target_density << ",\n"
        << "    \"compression_scale_per_dimension\": "
        << compression_scale << ",\n"
        << "    \"compression_axes\": \"x y z\",\n"
        << "    \"pressure_control\": \"isotropic xyz\",\n"
        << "    \"hot_temperature_K\": 800.0,\n"
        << "    \"final_temperature_K\": 300.0,\n"
        << "    \"timestep_fs\": " << kTimestepFs << ",\n"
        << "    \"equilibration_steps\": " << kEquilibrationSteps << ",\n"
        << "    \"bulk_300K_npt_steps\": " << kBulkFinalNptSteps << ",\n"
        << "    \"surface_tension_method\": \"film pressure anisotropy\",\n"
        << "    \"energy_sample_every_steps\": "
        << kEnergySampleEverySteps << ",\n"
        << "    \"total_run_steps\": " << kEquilibrationSteps << "\n"
        << "  },\n"
        << "  \"film_from_bulk\": {\n"
        << "    \"source_data\": \""
        << json_escape(files.bulk_equilibrated_data_basename) << "\",\n"
        << "    \"source_temperature_K\": 300.0,\n"
        << "    \"minimum_padding_each_z_face_angstrom\": " << padding << ",\n"
        << "    \"padding_source\": \""
        << (settings.film_padding > kMinimumFilmPaddingAngstrom
            ? "explicit" : "20 A safety minimum")
        << "\",\n"
        << "    \"initial_lz_rule\": \"bulk equilibrated Lz + 2 * actual padding; actual padding may increase to keep unwrapped chains clear of walls\",\n"
        << "    \"initial_lx_ly_rule\": \"inherit equilibrated bulk snapshot\",\n"
        << "    \"boundary_after_conversion\": \"p p f\",\n"
        << "    \"initial_wall_steps\": " << kFilmWallSteps << ",\n"
        << "    \"temporary_wall_steps\": " << kFilmWallSteps << ",\n"
        << "    \"wall_free_relaxation_steps\": 0,\n"
        << "    \"wall_free_relaxation_time_ns\": 0,\n"
        << "    \"guarded_relaxation_steps\": " << kFilmRelaxSteps << ",\n"
        << "    \"guarded_relaxation_time_ns\": "
        << kFilmRelaxSteps * kTimestepFs / 1.0e6 << ",\n"
        << "    \"film_production_steps\": "
        << kFilmProductionSteps << ",\n"
        << "    \"film_production_time_ns\": "
        << kFilmProductionSteps * kTimestepFs / 1.0e6 << ",\n"
        << "    \"total_run_steps\": "
        << kFilmWallSteps + kFilmRelaxSteps + kFilmProductionSteps << ",\n"
        << "    \"walls_during_production\": true,\n"
        << "    \"wall_force_output_columns\": [\"wall_lo_force\", \"wall_hi_force\"],\n"
        << "    \"valid_free_surface_requires_zero_wall_force\": true\n"
        << "  }\n"
        << "}\n";
    if (!out) throw std::runtime_error("Failed while writing info file: " + files.info);
}

void report(
    const Settings& settings,
    const OutputFiles& files,
    const System& system,
    const Box& box
) {
    const double realized_monomer_percent =
        100.0 * settings.total_mps_repeats / settings.total_repeats;
    const double realized_weight_percent =
        100.0 * settings.total_mps_repeats * kMpsRepeatMass /
        total_mass(settings);
    std::cerr << std::fixed << std::setprecision(4)
        << "Generated standalone silicone oil\n"
        << "  chains: " << settings.chains << ", repeat units/chain: ";
    if (settings.length_list_mode)
        std::cerr << "listed (" << settings.chain_count_entries.front().first
                  << '-' << settings.chain_count_entries.back().first << ")\n";
    else
        std::cerr << settings.length << '\n';
    std::cerr << "  MPS repeats/chain: ";
    if (settings.total_mps_repeats == 0)
        std::cerr << 0;
    else if (settings.length_list_mode)
        std::cerr << "allocated from the global composition";
    else
        std::cerr << settings.base_mps_per_chain;
    if (!settings.length_list_mode && settings.chains_with_extra_mps > 0)
        std::cerr << " or " << settings.base_mps_per_chain + 1
                  << " (" << settings.chains_with_extra_mps
                  << " chains with the extra MPS repeat)";
    std::cerr << '\n'
        << "  realized MPS: " << realized_monomer_percent
        << " monomer%, " << realized_weight_percent << " wt%\n"
        << "  atoms/bonds/angles/dihedrals: "
        << system.atoms.size() << '/' << system.bonds.size() << '/'
        << system.angles.size() << '/' << system.dihedrals.size() << '\n'
        << "  initial bulk box: " << box.lx << " x " << box.ly << " x "
        << box.lz << " A at "
        << settings.density << " g/cm^3\n"
        << "  whole chains inside primary box: yes (image flags 0 0 0)\n"
        << "  scripted 800 K compression target: "
        << settings.target_density << " g/cm^3\n"
        << "  bulk: " << kEquilibrationSteps * kTimestepFs / 1.0e6
        << " ns through 300 K NPT; no bulk production\n"
        << "  film: " << kFilmWallSteps * kTimestepFs / 1.0e6
        << " ns initialization, then " << kFilmRelaxSteps * kTimestepFs / 1.0e6
        << " ns equilibration and "
        << kFilmProductionSteps * kTimestepFs / 1.0e6
        << " ns pressure production with distant edge safety walls\n"
        << "  minimum film z padding: " << film_padding_angstrom(settings)
        << " A per face; check wall forces are zero during measurement\n"
        << "  generated files directory: " << files.directory << '\n'
        << "  submit both stages: bash " << files.pair_submit << '\n';
}

} // namespace

int main(int argc, char** argv) {
    try {
        Settings settings = parse_args(argc, argv);
        resolve_composition(settings);
        validate(settings);
        derive_output_name(settings);
        if (!settings.config_file.empty() &&
            std::filesystem::path(settings.output).is_relative()) {
            settings.output =
                (std::filesystem::path(settings.config_file).parent_path() /
                 settings.output).lexically_normal().string();
        }

        const Box box = calculate_box(settings);
        const System system = generate_system(settings, box);
        const OutputFiles files = output_files(settings);
        create_output_directory(files);
        write_data(files, system, box);
        write_input(settings, files);
        write_film_input(settings, files);
        write_submit(files);
        write_info(settings, files, system, box);
        report(settings, files, system, box);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
