#ifndef PDMS_NETWORK_COMMON_HPP
#define PDMS_NETWORK_COMMON_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pdms_analysis {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vec3 operator+(const Vec3 &a, const Vec3 &b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3 operator*(double scale, const Vec3 &value) {
    return {scale * value.x, scale * value.y, scale * value.z};
}
inline Vec3 &operator+=(Vec3 &a, const Vec3 &b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}
inline double norm2(const Vec3 &value) {
    return value.x * value.x + value.y * value.y + value.z * value.z;
}
inline double norm(const Vec3 &value) { return std::sqrt(norm2(value)); }

inline std::string trim(const std::string &input) {
    const std::size_t first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const std::size_t last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

inline bool begins_with(const std::string &text, const std::string &prefix) {
    return text.size() >= prefix.size() &&
           text.compare(0, prefix.size(), prefix) == 0;
}

inline std::string read_text_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open file: " + path);
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.good() && !input.eof())
        throw std::runtime_error("failed while reading file: " + path);
    return contents.str();
}

inline std::size_t matching_brace(
    const std::string &text, std::size_t open_position) {
    if (open_position >= text.size() || text[open_position] != '{')
        throw std::runtime_error("internal JSON object parsing error");
    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t i = open_position; i < text.size(); ++i) {
        const char c = text[i];
        if (in_string) {
            if (escaped) escaped = false;
            else if (c == '\\') escaped = true;
            else if (c == '"') in_string = false;
            continue;
        }
        if (c == '"') in_string = true;
        else if (c == '{') ++depth;
        else if (c == '}' && --depth == 0) return i;
    }
    throw std::runtime_error("unterminated JSON object");
}

inline std::string json_object_for_key(
    const std::string &text, const std::string &key) {
    const std::string quoted = "\"" + key + "\"";
    const std::size_t key_position = text.find(quoted);
    if (key_position == std::string::npos)
        throw std::runtime_error("missing JSON object: " + key);
    const std::size_t colon = text.find(':', key_position + quoted.size());
    const std::size_t open = text.find('{', colon + 1);
    if (colon == std::string::npos || open == std::string::npos)
        throw std::runtime_error("invalid JSON object: " + key);
    return text.substr(open, matching_brace(text, open) - open + 1);
}

inline std::string json_string(
    const std::string &text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern))
        throw std::runtime_error("missing string field in info file: " + key);
    return match[1].str();
}

inline double json_number(
    const std::string &text, const std::string &key) {
    const std::regex pattern(
        "\"" + key +
        "\"\\s*:\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern))
        throw std::runtime_error("missing numeric field in info file: " + key);
    return std::stod(match[1].str());
}

inline long long json_integer(
    const std::string &text, const std::string &key) {
    const double value = json_number(text, key);
    if (!std::isfinite(value) || std::fabs(value - std::round(value)) > 1.0e-8)
        throw std::runtime_error("info field is not an integer: " + key);
    return static_cast<long long>(std::llround(value));
}

inline std::vector<long long> json_integer_array(
    const std::string &text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]");
    std::smatch match;
    if (!std::regex_search(text, match, pattern))
        throw std::runtime_error("missing integer array in info file: " + key);
    std::vector<long long> result;
    std::string value;
    std::istringstream fields(match[1].str());
    while (std::getline(fields, value, ',')) {
        value = trim(value);
        if (!value.empty()) result.push_back(std::stoll(value));
    }
    return result;
}

struct ComponentInfo {
    std::string key;
    long long beads_per_molecule = 0;
    long long molecules = 0;
    long long molecule_start = 0;
    long long molecule_end = 0;
};

enum Component : int {
    kStrand = 0,
    kCrosslinker = 1,
    kModerator = 2,
    kFiller = 3,
    kComponentCount = 4
};

struct ModelInfo {
    std::string case_name;
    std::string geometry;
    std::string strand_topology;
    long long format_version = 0;
    double film_thickness_angstrom =
        std::numeric_limits<double>::quiet_NaN();
    double film_wall_cutoff_per_side_angstrom = 0.0;
    std::array<ComponentInfo, kComponentCount> components;
    long long total_beads = 0;
    long long total_molecules = 0;
    long long strand_functionality = 0;
    long long star_arm_count = 0;
    long long star_center_count = 0;
    long long star_arm_length = 0;
    long long graft_backbone_length = 0;
    long long graft_side_chain_length = 0;
    long long graft_spacing = 0;
    long long graft_side_chain_count = 0;
    std::vector<long long> reactive_bead_sites;
    std::vector<long long> crosslinker_reactive_bead_sites;
    int crosslink_bond_type = 2;
    double timestep_fs = 5.0;
    double final_temperature_k = 300.0;

    bool periodic_z() const { return geometry == "bulk"; }
};

inline ComponentInfo parse_component(
    const std::string &components, const std::string &key) {
    const std::string object = json_object_for_key(components, key);
    ComponentInfo result;
    result.key = key;
    result.beads_per_molecule = json_integer(object, "N");
    result.molecules = json_integer(object, "M");
    if (result.molecules > 0) {
        result.molecule_start = json_integer(object, "molecule_id_start");
        result.molecule_end = json_integer(object, "molecule_id_end");
    }
    return result;
}

inline long long optional_integer(
    const std::string &text, const std::string &key, long long fallback = 0) {
    const std::regex null_pattern("\"" + key + "\"\\s*:\\s*null");
    if (std::regex_search(text, null_pattern)) return fallback;
    try { return json_integer(text, key); }
    catch (const std::runtime_error &) { return fallback; }
}

inline double optional_number(
    const std::string &text, const std::string &key,
    double fallback = std::numeric_limits<double>::quiet_NaN()) {
    const std::regex null_pattern("\"" + key + "\"\\s*:\\s*null");
    if (std::regex_search(text, null_pattern)) return fallback;
    try { return json_number(text, key); }
    catch (const std::runtime_error &) { return fallback; }
}

inline ModelInfo parse_model_info(const std::string &path) {
    const std::string text = read_text_file(path);
    const std::string format = json_string(text, "format");
    ModelInfo info;
    info.format_version = json_integer(text, "format_version");
    if (format != "pdms-elastomer-model-info" || info.format_version < 3)
        throw std::runtime_error(
            "unsupported model info; expected pdms-elastomer-model-info version 3 or newer");
    info.case_name = json_string(text, "case_name");
    info.geometry = json_string(text, "geometry");
    if (info.geometry != "bulk" && info.geometry != "film")
        throw std::runtime_error("geometry must be bulk or film");
    if (info.geometry == "film") {
        info.film_thickness_angstrom = optional_number(
            text, "film_thickness_angstrom");
        info.film_wall_cutoff_per_side_angstrom = optional_number(
            text, "film_wall_cutoff_per_side_angstrom", 0.0);
        if (!(info.film_wall_cutoff_per_side_angstrom > 0.0)) {
            try {
                const std::string simulation =
                    json_object_for_key(text, "simulation_template");
                const std::string cold_lj =
                    json_object_for_key(simulation, "cold_lj");
                info.film_wall_cutoff_per_side_angstrom =
                    json_number(cold_lj, "cutoff");
            } catch (const std::runtime_error &) {
                info.film_wall_cutoff_per_side_angstrom = 0.0;
            }
        }
    }

    const std::string components = json_object_for_key(text, "components");
    info.components[kStrand] = parse_component(components, "strands");
    info.components[kCrosslinker] = parse_component(components, "crosslinkers");
    info.components[kModerator] = parse_component(components, "moderators");
    info.components[kFiller] = parse_component(components, "filler");
    const std::string composition = json_object_for_key(text, "composition");
    info.total_beads = json_integer(composition, "total_beads");
    info.total_molecules = json_integer(composition, "total_molecules");

    if (text.find("\"strand\"") != std::string::npos) {
        const std::string strand = json_object_for_key(text, "strand");
        info.strand_topology = json_string(strand, "topology");
        info.strand_functionality = json_integer(strand, "functionality");
        info.star_arm_count = optional_integer(strand, "star_arm_count");
        info.star_center_count = optional_integer(strand, "star_center_count");
        info.star_arm_length = optional_integer(strand, "star_arm_length");
        info.graft_backbone_length = optional_integer(strand, "grafted_backbone_length");
        info.graft_side_chain_length = optional_integer(strand, "grafted_side_chain_length");
        info.graft_spacing = optional_integer(strand, "graft_spacing");
        info.graft_side_chain_count = optional_integer(strand, "side_chain_count");
        info.reactive_bead_sites = json_integer_array(strand, "reactive_bead_sites");
    } else {
        // Version-3 linear files produced before architecture metadata was added.
        info.strand_topology = "linear";
        info.strand_functionality = 2;
        info.reactive_bead_sites = {
            1, info.components[kStrand].beads_per_molecule};
    }

    const std::string force_field = json_object_for_key(text, "force_field");
    const std::string bond_types = json_object_for_key(force_field, "bond_type_map");
    info.crosslink_bond_type = static_cast<int>(json_integer(bond_types, "crosslink"));
    const std::string crosslinker = json_object_for_key(text, "crosslinker");
    info.crosslinker_reactive_bead_sites =
        json_integer_array(crosslinker, "reactive_bead_sites");
    const std::string simulation = json_object_for_key(text, "simulation_template");
    info.timestep_fs = json_number(simulation, "timestep_fs");
    info.final_temperature_k = json_number(simulation, "final_temperature_K");
    return info;
}

inline int component_for_molecule(long long molecule, const ModelInfo &info) {
    for (int component = 0; component < kComponentCount; ++component) {
        const ComponentInfo &entry = info.components[static_cast<std::size_t>(component)];
        if (entry.molecules > 0 && molecule >= entry.molecule_start &&
            molecule <= entry.molecule_end)
            return component;
    }
    throw std::runtime_error(
        "molecule ID outside component ranges: " + std::to_string(molecule));
}

struct Box {
    double xlo = 0.0, xhi = 0.0;
    double ylo = 0.0, yhi = 0.0;
    double zlo = 0.0, zhi = 0.0;
    bool have_x = false, have_y = false, have_z = false;
    double lx() const { return xhi - xlo; }
    double ly() const { return yhi - ylo; }
    double lz() const { return zhi - zlo; }
};

struct Atom {
    bool seen = false;
    long long id = 0;
    long long molecule = 0;
    int type = 0;
    double mass = 0.0;
    Vec3 position;
    Vec3 velocity;
    bool have_velocity = false;
    std::array<long long, 3> image{{0, 0, 0}};
    bool have_images = false;
};

struct Bond {
    long long id = 0;
    int type = 0;
    long long first = 0;
    long long second = 0;
};

struct DataFile {
    Box box;
    long long declared_atoms = 0;
    long long declared_bonds = 0;
    long long declared_angles = 0;
    long long declared_dihedrals = 0;
    long long velocities_read = 0;
    double total_mass_g_per_mol = 0.0;
    std::vector<double> masses;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<std::vector<long long>> molecule_atoms;
};

enum class DataSection {
    kHeader, kMasses, kAtoms, kVelocities, kBonds, kOther
};

inline bool section_line(const std::string &line, const std::string &name) {
    return begins_with(line, name) &&
        (line.size() == name.size() ||
         std::isspace(static_cast<unsigned char>(line[name.size()])) ||
         line[name.size()] == '#');
}

inline void parse_header_line(const std::string &line, DataFile &data) {
    long long count = 0;
    std::string first, second;
    std::istringstream counts(line);
    if (counts >> count >> first) {
        if (first == "atoms") data.declared_atoms = count;
        else if (first == "bonds") data.declared_bonds = count;
        else if (first == "angles") data.declared_angles = count;
        else if (first == "dihedrals") data.declared_dihedrals = count;
    }
    double lo = 0.0, hi = 0.0;
    std::string lower, upper;
    std::istringstream bounds(line);
    if (bounds >> lo >> hi >> lower >> upper) {
        if (lower == "xlo" && upper == "xhi") {
            data.box.xlo = lo; data.box.xhi = hi; data.box.have_x = true;
        } else if (lower == "ylo" && upper == "yhi") {
            data.box.ylo = lo; data.box.yhi = hi; data.box.have_y = true;
        } else if (lower == "zlo" && upper == "zhi") {
            data.box.zlo = lo; data.box.zhi = hi; data.box.have_z = true;
        }
    }
}

inline DataFile parse_data_file(
    const std::string &path, const ModelInfo &info) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open LAMMPS data file: " + path);
    DataFile data;
    data.molecule_atoms.resize(static_cast<std::size_t>(info.total_molecules + 1));
    DataSection section = DataSection::kHeader;
    std::string line;
    long long line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string clean = trim(line);
        if (clean.empty() || clean.front() == '#') continue;
        if (section_line(clean, "Masses")) {
            section = DataSection::kMasses;
            continue;
        }
        if (section_line(clean, "Atoms")) {
            if (data.declared_atoms <= 0)
                throw std::runtime_error("atom count must precede Atoms section");
            data.atoms.resize(static_cast<std::size_t>(data.declared_atoms + 1));
            section = DataSection::kAtoms;
            continue;
        }
        if (section_line(clean, "Bonds")) {
            section = DataSection::kBonds;
            continue;
        }
        if (section_line(clean, "Velocities")) {
            section = DataSection::kVelocities;
            continue;
        }
        if (section_line(clean, "Angles") || section_line(clean, "Dihedrals") ||
            section_line(clean, "Pair Coeffs") || section_line(clean, "Bond Coeffs") ||
            section_line(clean, "Angle Coeffs") || section_line(clean, "Dihedral Coeffs") ||
            section_line(clean, "Impropers")) {
            section = DataSection::kOther;
            continue;
        }
        if (section == DataSection::kHeader) {
            parse_header_line(clean, data);
        } else if (section == DataSection::kMasses) {
            int type = 0;
            double mass = 0.0;
            std::istringstream fields(clean);
            if (!(fields >> type >> mass) || type < 1 || mass <= 0.0)
                throw std::runtime_error(
                    "invalid mass at line " + std::to_string(line_number));
            if (data.masses.size() <= static_cast<std::size_t>(type))
                data.masses.resize(static_cast<std::size_t>(type + 1), 0.0);
            if (data.masses[static_cast<std::size_t>(type)] != 0.0)
                throw std::runtime_error("duplicate atom-type mass");
            data.masses[static_cast<std::size_t>(type)] = mass;
        } else if (section == DataSection::kAtoms) {
            Atom atom;
            double charge = 0.0;
            std::istringstream fields(clean);
            if (!(fields >> atom.id >> atom.molecule >> atom.type >> charge >>
                  atom.position.x >> atom.position.y >> atom.position.z))
                throw std::runtime_error("invalid atom at line " + std::to_string(line_number));
            atom.have_images = static_cast<bool>(
                fields >> atom.image[0] >> atom.image[1] >> atom.image[2]);
            if (atom.id < 1 || atom.id > data.declared_atoms ||
                data.atoms[static_cast<std::size_t>(atom.id)].seen)
                throw std::runtime_error("invalid or duplicate atom ID");
            component_for_molecule(atom.molecule, info);
            atom.seen = true;
            data.atoms[static_cast<std::size_t>(atom.id)] = atom;
            data.molecule_atoms[static_cast<std::size_t>(atom.molecule)].push_back(atom.id);
        } else if (section == DataSection::kVelocities) {
            long long id = 0;
            Vec3 velocity;
            std::istringstream fields(clean);
            if (!(fields >> id >> velocity.x >> velocity.y >> velocity.z) ||
                id < 1 || id > data.declared_atoms || data.atoms.empty() ||
                !data.atoms[static_cast<std::size_t>(id)].seen)
                throw std::runtime_error(
                    "invalid velocity at line " + std::to_string(line_number));
            Atom &atom = data.atoms[static_cast<std::size_t>(id)];
            if (atom.have_velocity)
                throw std::runtime_error("duplicate atom velocity");
            atom.velocity = velocity;
            atom.have_velocity = true;
            ++data.velocities_read;
        } else if (section == DataSection::kBonds) {
            Bond bond;
            std::istringstream fields(clean);
            if (!(fields >> bond.id >> bond.type >> bond.first >> bond.second))
                throw std::runtime_error("invalid bond at line " + std::to_string(line_number));
            if (bond.first < 1 || bond.first > data.declared_atoms ||
                bond.second < 1 || bond.second > data.declared_atoms)
                throw std::runtime_error("bond atom ID outside atom range");
            data.bonds.push_back(bond);
        }
    }
    if (!data.box.have_x || !data.box.have_y || !data.box.have_z)
        throw std::runtime_error("incomplete orthorhombic box bounds");
    if (static_cast<long long>(data.bonds.size()) != data.declared_bonds)
        throw std::runtime_error("bond count does not match LAMMPS header");
    for (long long id = 1; id <= data.declared_atoms; ++id)
        if (!data.atoms[static_cast<std::size_t>(id)].seen)
            throw std::runtime_error("missing atom ID " + std::to_string(id));
        else {
            Atom &atom = data.atoms[static_cast<std::size_t>(id)];
            if (atom.type < 1 ||
                static_cast<std::size_t>(atom.type) >= data.masses.size() ||
                data.masses[static_cast<std::size_t>(atom.type)] <= 0.0)
                throw std::runtime_error("missing mass for atom type " +
                                         std::to_string(atom.type));
            atom.mass = data.masses[static_cast<std::size_t>(atom.type)];
            data.total_mass_g_per_mol += atom.mass;
        }
    for (auto &atoms : data.molecule_atoms) std::sort(atoms.begin(), atoms.end());
    return data;
}

inline double minimum_image(double delta, double length) {
    return delta - std::round(delta / length) * length;
}

inline Vec3 minimum_image_vector(
    const Vec3 &delta, const Box &box, bool periodic_z) {
    return {
        minimum_image(delta.x, box.lx()),
        minimum_image(delta.y, box.ly()),
        periodic_z ? minimum_image(delta.z, box.lz()) : delta.z
    };
}

inline Vec3 wrap_position(Vec3 value, const Box &box, bool periodic_z) {
    const auto wrap = [](double x, double lo, double length) {
        x -= std::floor((x - lo) / length) * length;
        if (x >= lo + length) x -= length;
        return x;
    };
    value.x = wrap(value.x, box.xlo, box.lx());
    value.y = wrap(value.y, box.ylo, box.ly());
    if (periodic_z) value.z = wrap(value.z, box.zlo, box.lz());
    return value;
}

class DisjointSet {
  public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }
    long long find(long long value) {
        long long &parent = parent_.at(static_cast<std::size_t>(value));
        if (parent != value) parent = find(parent);
        return parent;
    }
    void unite(long long first, long long second) {
        first = find(first); second = find(second);
        if (first == second) return;
        if (rank_[static_cast<std::size_t>(first)] <
            rank_[static_cast<std::size_t>(second)]) std::swap(first, second);
        parent_[static_cast<std::size_t>(second)] = first;
        if (rank_[static_cast<std::size_t>(first)] ==
            rank_[static_cast<std::size_t>(second)])
            ++rank_[static_cast<std::size_t>(first)];
    }
  private:
    std::vector<long long> parent_;
    std::vector<unsigned char> rank_;
};

struct JunctionNode {
    long long id = 0;
    std::string kind;
    long long representative_molecule = 0;
    std::vector<long long> member_molecules;
    std::vector<long long> attachment_atoms;
    Vec3 position;
    long long active_degree = 0;
    long long dangling_degree = 0;
    long long component = 0;
};

struct EffectiveStrand {
    long long id = 0;
    long long parent_molecule = 0;
    std::string parent_topology;
    std::string status;
    long long first_node = 0;
    long long second_node = 0;
    long long first_atom = 0;
    long long second_atom = 0;
    std::vector<long long> atoms;
    long long contour_bonds = 0;
    double contour_length = 0.0;
    Vec3 end_to_end;
    double end_to_end_length = 0.0;
    double straightness = 0.0;
    std::array<int, 3> winding{{0, 0, 0}};
    long long graph_component = 0;
};

struct ParentRecord {
    long long molecule = 0;
    std::string topology;
    long long possible_sites = 0;
    long long reacted_sites = 0;
    long long effective_strands = 0;
    std::string state;
};

struct ReducedNetwork {
    std::vector<JunctionNode> nodes{{}};
    std::vector<EffectiveStrand> strands;
    std::vector<ParentRecord> parents;
    std::unordered_map<long long, long long> network_atom_node;
    long long strand_reactive_total = 0;
    long long strand_reactive_reacted = 0;
    long long crosslinker_reactive_total = 0;
    long long crosslinker_reactive_reacted = 0;
};

inline std::vector<std::vector<long long>> strand_internal_adjacency(
    const DataFile &data, const ModelInfo &info) {
    std::vector<std::vector<long long>> adjacency(data.atoms.size());
    for (const Bond &bond : data.bonds) {
        if (bond.type == info.crosslink_bond_type) continue;
        const Atom &first = data.atoms[static_cast<std::size_t>(bond.first)];
        const Atom &second = data.atoms[static_cast<std::size_t>(bond.second)];
        if (first.molecule != second.molecule) continue;
        if (component_for_molecule(first.molecule, info) != kStrand) continue;
        adjacency[static_cast<std::size_t>(bond.first)].push_back(bond.second);
        adjacency[static_cast<std::size_t>(bond.second)].push_back(bond.first);
    }
    for (auto &neighbors : adjacency) std::sort(neighbors.begin(), neighbors.end());
    return adjacency;
}

inline std::vector<long long> unique_path(
    long long start, long long finish,
    const std::vector<std::vector<long long>> &adjacency,
    const std::set<std::pair<long long, long long>> &blocked = {}) {
    std::vector<long long> parent(adjacency.size(), -1);
    std::queue<long long> queue;
    parent[static_cast<std::size_t>(start)] = start;
    queue.push(start);
    while (!queue.empty() && parent[static_cast<std::size_t>(finish)] < 0) {
        const long long current = queue.front(); queue.pop();
        for (long long neighbor : adjacency[static_cast<std::size_t>(current)]) {
            if (blocked.count(std::minmax(current, neighbor))) continue;
            if (parent[static_cast<std::size_t>(neighbor)] >= 0) continue;
            parent[static_cast<std::size_t>(neighbor)] = current;
            queue.push(neighbor);
        }
    }
    if (parent[static_cast<std::size_t>(finish)] < 0)
        throw std::runtime_error("no intramolecular path between selected atoms");
    std::vector<long long> path;
    for (long long atom = finish;; atom = parent[static_cast<std::size_t>(atom)]) {
        path.push_back(atom);
        if (atom == start) break;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

inline std::vector<long long> linear_order(
    const std::vector<long long> &atoms,
    const std::vector<std::vector<long long>> &adjacency) {
    std::vector<long long> ends;
    for (long long atom : atoms)
        if (adjacency[static_cast<std::size_t>(atom)].size() == 1) ends.push_back(atom);
    if (ends.size() != 2)
        throw std::runtime_error("linear strand does not have exactly two ends");
    return unique_path(ends[0], ends[1], adjacency);
}

inline std::vector<long long> ring_order(
    const std::vector<long long> &atoms,
    const std::vector<std::vector<long long>> &adjacency) {
    if (atoms.size() < 3) throw std::runtime_error("ring requires at least three beads");
    for (long long atom : atoms)
        if (adjacency[static_cast<std::size_t>(atom)].size() != 2)
            throw std::runtime_error("ring molecule is not a simple cycle");
    const long long start = atoms.front();
    long long previous = -1;
    long long current = start;
    std::vector<long long> order;
    do {
        order.push_back(current);
        const auto &neighbors = adjacency[static_cast<std::size_t>(current)];
        const long long next = neighbors[0] == previous ? neighbors[1] : neighbors[0];
        previous = current;
        current = next;
        if (order.size() > atoms.size())
            throw std::runtime_error("ring traversal did not close");
    } while (current != start);
    if (order.size() != atoms.size())
        throw std::runtime_error("ring traversal omitted beads");
    return order;
}

inline std::vector<Vec3> unwrapped_path(
    const std::vector<long long> &path, const DataFile &data,
    const ModelInfo &info) {
    std::vector<Vec3> positions;
    if (path.empty()) return positions;
    positions.push_back(data.atoms[static_cast<std::size_t>(path.front())].position);
    for (std::size_t i = 1; i < path.size(); ++i) {
        const Vec3 previous = data.atoms[static_cast<std::size_t>(path[i - 1])].position;
        const Vec3 current = data.atoms[static_cast<std::size_t>(path[i])].position;
        positions.push_back(positions.back() + minimum_image_vector(
            current - previous, data.box, info.periodic_z()));
    }
    return positions;
}

inline std::string classify_edge(long long first_node, long long second_node) {
    if (first_node == 0 && second_node == 0) return "isolated";
    if (first_node == 0 || second_node == 0) return "dangling";
    if (first_node == second_node) return "self_loop";
    return "active";
}

inline EffectiveStrand make_effective_strand(
    long long parent, const std::string &topology,
    const std::vector<long long> &path,
    long long first_node, long long second_node,
    const DataFile &data, const ModelInfo &info,
    const std::string &forced_status = "") {
    if (path.size() < 2)
        throw std::runtime_error("effective strand contains fewer than two contour points");
    EffectiveStrand strand;
    strand.parent_molecule = parent;
    strand.parent_topology = topology;
    strand.first_node = first_node;
    strand.second_node = second_node;
    strand.first_atom = path.front();
    strand.second_atom = path.back();
    strand.atoms = path;
    strand.contour_bonds = static_cast<long long>(path.size() - 1);
    strand.status = forced_status.empty()
        ? classify_edge(first_node, second_node) : forced_status;
    const std::vector<Vec3> positions = unwrapped_path(path, data, info);
    for (std::size_t i = 1; i < positions.size(); ++i)
        strand.contour_length += norm(positions[i] - positions[i - 1]);
    strand.end_to_end = positions.back() - positions.front();
    strand.end_to_end_length = norm(strand.end_to_end);
    if (strand.contour_length > 0.0)
        strand.straightness = strand.end_to_end_length / strand.contour_length;
    return strand;
}

inline long long local_rank(
    long long atom, const std::vector<long long> &molecule_atoms) {
    const auto found = std::lower_bound(molecule_atoms.begin(), molecule_atoms.end(), atom);
    if (found == molecule_atoms.end() || *found != atom)
        throw std::runtime_error("atom missing from its molecule list");
    return static_cast<long long>(found - molecule_atoms.begin()) + 1;
}

inline ReducedNetwork reduce_network(
    const DataFile &data, const ModelInfo &info) {
    ReducedNetwork reduced;
    const auto internal = strand_internal_adjacency(data, info);
    DisjointSet junction_sets(static_cast<std::size_t>(info.total_molecules + 1));

    // Moderator-to-crosslinker reactions make one collapsed junction cluster.
    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        const long long first_molecule =
            data.atoms[static_cast<std::size_t>(bond.first)].molecule;
        const long long second_molecule =
            data.atoms[static_cast<std::size_t>(bond.second)].molecule;
        const int first_component = component_for_molecule(first_molecule, info);
        const int second_component = component_for_molecule(second_molecule, info);
        const bool first_junction =
            first_component == kCrosslinker || first_component == kModerator;
        const bool second_junction =
            second_component == kCrosslinker || second_component == kModerator;
        if (first_junction && second_junction)
            junction_sets.unite(first_molecule, second_molecule);
    }

    std::map<long long, std::vector<long long>> junction_clusters;
    for (int component : {kCrosslinker, kModerator}) {
        const ComponentInfo &entry = info.components[static_cast<std::size_t>(component)];
        for (long long molecule = entry.molecule_start;
             molecule > 0 && molecule <= entry.molecule_end; ++molecule)
            junction_clusters[junction_sets.find(molecule)].push_back(molecule);
    }
    std::unordered_map<long long, long long> molecule_node;
    auto ensure_chemical_node = [&](long long molecule) {
        const auto existing = molecule_node.find(molecule);
        if (existing != molecule_node.end()) return existing->second;
        const long long root = junction_sets.find(molecule);
        const auto cluster = junction_clusters.find(root);
        if (cluster == junction_clusters.end())
            throw std::runtime_error("reaction partner is not a junction molecule");
        JunctionNode node;
        node.id = static_cast<long long>(reduced.nodes.size());
        node.kind = "chemical_junction";
        node.representative_molecule = cluster->second.front();
        node.member_molecules = cluster->second;
        reduced.nodes.push_back(node);
        for (long long member : cluster->second) molecule_node[member] = node.id;
        return node.id;
    };

    // Map actually reacted component-1 sites to their collapsed chemical junction.
    std::set<long long> reacted_crosslinker_atoms;
    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        const Atom &first = data.atoms[static_cast<std::size_t>(bond.first)];
        const Atom &second = data.atoms[static_cast<std::size_t>(bond.second)];
        const int first_component = component_for_molecule(first.molecule, info);
        const int second_component = component_for_molecule(second.molecule, info);
        if (first_component == kCrosslinker) reacted_crosslinker_atoms.insert(first.id);
        if (second_component == kCrosslinker) reacted_crosslinker_atoms.insert(second.id);
        const Atom *network_atom = nullptr;
        const Atom *junction_atom = nullptr;
        if (first_component == kStrand &&
            (second_component == kCrosslinker || second_component == kModerator)) {
            network_atom = &first; junction_atom = &second;
        } else if (second_component == kStrand &&
                   (first_component == kCrosslinker || first_component == kModerator)) {
            network_atom = &second; junction_atom = &first;
        }
        if (network_atom == nullptr) continue;
        const long long node_id = ensure_chemical_node(junction_atom->molecule);
        const auto inserted = reduced.network_atom_node.emplace(network_atom->id, node_id);
        if (!inserted.second && inserted.first->second != node_id)
            throw std::runtime_error(
                "one strand reactive bead is bonded to multiple junction clusters");
        reduced.nodes[static_cast<std::size_t>(node_id)]
            .attachment_atoms.push_back(network_atom->id);
    }

    const ComponentInfo &strand_component = info.components[kStrand];
    for (long long molecule = strand_component.molecule_start;
         molecule > 0 && molecule <= strand_component.molecule_end; ++molecule) {
        const auto &atoms = data.molecule_atoms[static_cast<std::size_t>(molecule)];
        if (atoms.empty()) throw std::runtime_error("empty strand molecule");
        std::vector<long long> possible_sites;
        std::vector<long long> reacted_sites;
        const std::set<long long> declared_sites(
            info.reactive_bead_sites.begin(), info.reactive_bead_sites.end());
        for (long long atom : atoms) {
            if (!declared_sites.count(local_rank(atom, atoms))) continue;
            possible_sites.push_back(atom);
            if (reduced.network_atom_node.count(atom)) reacted_sites.push_back(atom);
        }
        reduced.strand_reactive_total += static_cast<long long>(possible_sites.size());
        reduced.strand_reactive_reacted += static_cast<long long>(reacted_sites.size());
        ParentRecord parent;
        parent.molecule = molecule;
        parent.topology = info.strand_topology;
        parent.possible_sites = static_cast<long long>(possible_sites.size());
        parent.reacted_sites = static_cast<long long>(reacted_sites.size());
        const std::size_t before = reduced.strands.size();

        if (info.strand_topology == "linear") {
            const std::vector<long long> path = linear_order(atoms, internal);
            const long long first_node = reduced.network_atom_node.count(path.front())
                ? reduced.network_atom_node.at(path.front()) : 0;
            const long long second_node = reduced.network_atom_node.count(path.back())
                ? reduced.network_atom_node.at(path.back()) : 0;
            reduced.strands.push_back(make_effective_strand(
                molecule, "linear", path, first_node, second_node, data, info));
            parent.state = classify_edge(first_node, second_node);
        } else if (info.strand_topology == "ring") {
            const std::vector<long long> cycle = ring_order(atoms, internal);
            std::map<long long, std::size_t> cycle_position;
            for (std::size_t i = 0; i < cycle.size(); ++i) cycle_position[cycle[i]] = i;
            std::sort(reacted_sites.begin(), reacted_sites.end(),
                [&](long long a, long long b) {
                    return cycle_position.at(a) < cycle_position.at(b);
                });
            if (reacted_sites.empty()) {
                parent.state = "isolated_ring";
            } else if (reacted_sites.size() == 1) {
                const std::size_t start = cycle_position.at(reacted_sites.front());
                std::vector<long long> path;
                for (std::size_t step = 0; step <= cycle.size(); ++step)
                    path.push_back(cycle[(start + step) % cycle.size()]);
                const long long node = reduced.network_atom_node.at(reacted_sites.front());
                reduced.strands.push_back(make_effective_strand(
                    molecule, "ring", path, node, node, data, info,
                    "dangling_loop"));
                parent.state = "dangling_ring";
            } else {
                for (std::size_t site = 0; site < reacted_sites.size(); ++site) {
                    const std::size_t start = cycle_position.at(reacted_sites[site]);
                    const std::size_t finish = cycle_position.at(
                        reacted_sites[(site + 1) % reacted_sites.size()]);
                    std::vector<long long> path;
                    for (std::size_t position = start;;
                         position = (position + 1) % cycle.size()) {
                        path.push_back(cycle[position]);
                        if (position == finish) break;
                    }
                    reduced.strands.push_back(make_effective_strand(
                        molecule, "ring", path,
                        reduced.network_atom_node.at(path.front()),
                        reduced.network_atom_node.at(path.back()), data, info));
                }
                parent.state = reacted_sites.size() == possible_sites.size()
                    ? "fully_reacted_ring" : "partially_reacted_ring";
            }
        } else if (info.strand_topology == "star") {
            if (info.star_center_count < 1 ||
                info.star_center_count >= static_cast<long long>(atoms.size()))
                throw std::runtime_error("invalid star center count in info file");
            JunctionNode center;
            center.id = static_cast<long long>(reduced.nodes.size());
            center.kind = "star_center";
            center.representative_molecule = molecule;
            center.member_molecules.push_back(molecule);
            std::vector<long long> centers;
            for (long long index = 0; index < info.star_center_count; ++index) {
                centers.push_back(atoms[static_cast<std::size_t>(index)]);
                center.attachment_atoms.push_back(atoms[static_cast<std::size_t>(index)]);
            }
            reduced.nodes.push_back(center);
            for (long long site : possible_sites) {
                std::vector<long long> best_path;
                for (long long core : centers) {
                    const std::vector<long long> candidate = unique_path(core, site, internal);
                    if (best_path.empty() || candidate.size() < best_path.size())
                        best_path = candidate;
                }
                const long long end_node = reduced.network_atom_node.count(site)
                    ? reduced.network_atom_node.at(site) : 0;
                reduced.strands.push_back(make_effective_strand(
                    molecule, "star", best_path, center.id, end_node, data, info));
            }
            parent.state = reacted_sites.empty() ? "isolated_star" :
                (reacted_sites.size() == possible_sites.size()
                    ? "fully_reacted_star" : "partially_reacted_star");
        } else if (info.strand_topology == "grafted") {
            std::sort(reacted_sites.begin(), reacted_sites.end(),
                [&](long long a, long long b) {
                    return local_rank(a, atoms) < local_rank(b, atoms);
                });
            if (reacted_sites.size() < 2) {
                parent.state = reacted_sites.empty()
                    ? "isolated_grafted" : "single_anchor_grafted";
            } else {
                for (std::size_t site = 1; site < reacted_sites.size(); ++site) {
                    const long long first = reacted_sites[site - 1];
                    const long long second = reacted_sites[site];
                    reduced.strands.push_back(make_effective_strand(
                        molecule, "grafted", unique_path(first, second, internal),
                        reduced.network_atom_node.at(first),
                        reduced.network_atom_node.at(second), data, info));
                }
                parent.state = reacted_sites.size() == possible_sites.size()
                    ? "fully_reacted_grafted" : "partially_reacted_grafted";
            }
        } else {
            throw std::runtime_error(
                "unsupported strand topology: " + info.strand_topology);
        }
        parent.effective_strands =
            static_cast<long long>(reduced.strands.size() - before);
        reduced.parents.push_back(parent);
    }

    const ComponentInfo &crosslinkers = info.components[kCrosslinker];
    const std::set<long long> declared_crosslinker_sites(
        info.crosslinker_reactive_bead_sites.begin(),
        info.crosslinker_reactive_bead_sites.end());
    for (long long molecule = crosslinkers.molecule_start;
         molecule > 0 && molecule <= crosslinkers.molecule_end; ++molecule) {
        const auto &atoms = data.molecule_atoms[static_cast<std::size_t>(molecule)];
        for (long long atom : atoms) {
            if (!declared_crosslinker_sites.count(local_rank(atom, atoms))) continue;
            ++reduced.crosslinker_reactive_total;
            if (reacted_crosslinker_atoms.count(atom))
                ++reduced.crosslinker_reactive_reacted;
        }
    }

    for (EffectiveStrand &strand : reduced.strands)
        strand.id = &strand - reduced.strands.data() + 1;

    // Junction positions use a locally unwrapped mean of their attachment points.
    for (std::size_t index = 1; index < reduced.nodes.size(); ++index) {
        JunctionNode &node = reduced.nodes[index];
        if (node.attachment_atoms.empty()) continue;
        const Vec3 reference = data.atoms[
            static_cast<std::size_t>(node.attachment_atoms.front())].position;
        Vec3 sum = reference;
        for (std::size_t i = 1; i < node.attachment_atoms.size(); ++i) {
            const Vec3 position = data.atoms[
                static_cast<std::size_t>(node.attachment_atoms[i])].position;
            sum += reference + minimum_image_vector(
                position - reference, data.box, info.periodic_z());
        }
        node.position = wrap_position(
            (1.0 / node.attachment_atoms.size()) * sum,
            data.box, info.periodic_z());
    }

    // Assign winding shifts after node positions are available.
    for (EffectiveStrand &strand : reduced.strands) {
        if (strand.first_node == 0 || strand.second_node == 0 ||
            strand.first_node == strand.second_node) continue;
        const JunctionNode &first = reduced.nodes[
            static_cast<std::size_t>(strand.first_node)];
        const JunctionNode &second = reduced.nodes[
            static_cast<std::size_t>(strand.second_node)];
        const Vec3 first_site = data.atoms[
            static_cast<std::size_t>(strand.first_atom)].position;
        const Vec3 second_site = data.atoms[
            static_cast<std::size_t>(strand.second_atom)].position;
        const Vec3 first_offset = minimum_image_vector(
            first_site - first.position, data.box, info.periodic_z());
        const Vec3 second_offset = minimum_image_vector(
            second_site - second.position, data.box, info.periodic_z());
        const Vec3 center_displacement =
            first_offset + strand.end_to_end - second_offset;
        const Vec3 base = minimum_image_vector(
            second.position - first.position, data.box, info.periodic_z());
        strand.winding[0] = static_cast<int>(
            std::llround((center_displacement.x - base.x) / data.box.lx()));
        strand.winding[1] = static_cast<int>(
            std::llround((center_displacement.y - base.y) / data.box.ly()));
        if (info.periodic_z())
            strand.winding[2] = static_cast<int>(
                std::llround((center_displacement.z - base.z) / data.box.lz()));
    }

    // Active graph connected components and degrees.
    DisjointSet graph_sets(reduced.nodes.size());
    for (const EffectiveStrand &strand : reduced.strands) {
        if (strand.status == "active") {
            graph_sets.unite(strand.first_node, strand.second_node);
            ++reduced.nodes[static_cast<std::size_t>(strand.first_node)].active_degree;
            ++reduced.nodes[static_cast<std::size_t>(strand.second_node)].active_degree;
        } else if (strand.status == "self_loop") {
            reduced.nodes[static_cast<std::size_t>(strand.first_node)].active_degree += 2;
        } else if (strand.status == "dangling" || strand.status == "dangling_loop") {
            const long long node = strand.first_node != 0
                ? strand.first_node : strand.second_node;
            if (node != 0)
                ++reduced.nodes[static_cast<std::size_t>(node)].dangling_degree;
        }
    }
    std::map<long long, long long> component_labels;
    for (std::size_t node = 1; node < reduced.nodes.size(); ++node) {
        if (reduced.nodes[node].active_degree == 0) continue;
        const long long root = graph_sets.find(static_cast<long long>(node));
        auto inserted = component_labels.emplace(
            root, static_cast<long long>(component_labels.size() + 1));
        reduced.nodes[node].component = inserted.first->second;
    }
    for (EffectiveStrand &strand : reduced.strands) {
        if (strand.first_node != 0)
            strand.graph_component = reduced.nodes[
                static_cast<std::size_t>(strand.first_node)].component;
        else if (strand.second_node != 0)
            strand.graph_component = reduced.nodes[
                static_cast<std::size_t>(strand.second_node)].component;
    }
    return reduced;
}

inline std::filesystem::path analysis_directory(
    const std::string &data_path, const ModelInfo &info,
    const std::string &requested = "") {
    if (!requested.empty()) return std::filesystem::path(requested);
    return std::filesystem::path(data_path).parent_path() /
           ("analysis_" + info.case_name);
}

inline void create_directory(const std::filesystem::path &path) {
    std::error_code error;
    std::filesystem::create_directories(path, error);
    if (error) throw std::runtime_error(
        "cannot create output directory " + path.string() + ": " + error.message());
}

inline std::string safe_case_name(std::string name) {
    for (char &c : name)
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '-' && c != '_') c = '_';
    return name.empty() ? "case" : name;
}

} // namespace pdms_analysis

#endif
