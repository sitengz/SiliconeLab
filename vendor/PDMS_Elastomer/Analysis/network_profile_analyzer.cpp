#include "network_common.hpp"

#include <iostream>

namespace {
using namespace pdms_analysis;

struct Options {
    std::string data_file;
    std::string info_file;
    std::string debye_waller_trajectory_file;
    std::string trajectory_file;
    std::string z1_sp_file;
    bool disable_z1 = false;
    std::string output_directory;
    double bin_width = 5.0;
    double profile_window_width = 5.0;
    double profile_step = 1.0;
    double wall_fraction = 0.20;
    double core_fraction = 0.20;
    double diffusion_fit_start_fraction = 0.50;
    double debye_waller_time_ps = std::numeric_limits<double>::quiet_NaN();
    double debye_waller_search_start_ps = 0.50;
    double debye_waller_search_end_ps = 20.0;
    long long frame_stride = 1;
    bool time_averaged_msd = false;
    long long time_origin_stride = 100;
    long long time_origin_count = 10;
    bool time_averaged_dw = false;
    long long dw_time_origin_stride = 100;
    long long dw_time_origin_count = 10;
};

struct SlidingGrid {
    int fine_bins = 1;
    double step_A = 1.0;
    int window_bins = 1;
    int half_window_bins = 0;
    double window_width_A = 1.0;
};

SlidingGrid make_sliding_grid(const Box &box, const Options &options) {
    SlidingGrid grid;
    grid.fine_bins = std::max(1, static_cast<int>(
        std::ceil(box.lz() / options.profile_step)));
    grid.step_A = box.lz() / grid.fine_bins;
    const double target_window_bins =
        options.profile_window_width / grid.step_A;
    const int maximum_odd_bins =
        grid.fine_bins % 2 == 0 ? grid.fine_bins - 1 : grid.fine_bins;
    grid.window_bins = std::clamp(
        static_cast<int>(std::llround(target_window_bins)), 1,
        maximum_odd_bins);
    if (grid.window_bins % 2 == 0) {
        const int lower = grid.window_bins - 1;
        const int upper = grid.window_bins + 1;
        if (upper <= maximum_odd_bins &&
            std::abs(upper - target_window_bins) <
                std::abs(lower - target_window_bins))
            grid.window_bins = upper;
        else
            grid.window_bins = lower;
    }
    grid.half_window_bins = grid.window_bins / 2;
    grid.window_width_A = grid.window_bins * grid.step_A;
    return grid;
}

constexpr double kAvogadroAngstrom = 0.602214076;

struct BinCounts {
    long long crosslink_bonds = 0;
    long long strand_crosslink_bonds = 0;
    long long moderator_crosslink_bonds = 0;
    long long junctions = 0;
    long long active_strands = 0;
    long long dangling_ends = 0;
    long long dangling_loops = 0;
    long long self_loops = 0;
    long long isolated_parents = 0;
};

struct BinProfile {
    BinCounts topology;
    std::array<long long, kComponentCount> beads{{0, 0, 0, 0}};
    std::array<double, kComponentCount> mass{{0.0, 0.0, 0.0, 0.0}};
    // Strand, crosslinker, and moderator functional sites, respectively.
    std::array<long long, 3> functional_total{{0, 0, 0}};
    std::array<long long, 3> functional_reacted{{0, 0, 0}};

    long long strand_markers = 0;
    long long strand_orientation_count = 0;
    double strand_ree_x2_sum = 0.0;
    double strand_ree_y2_sum = 0.0;
    double strand_ree_z2_sum = 0.0;
    double strand_stretch_sum = 0.0;
    double strand_p2z_sum = 0.0;
    long long active_orientation_count = 0;
    double active_ree_x2_sum = 0.0;
    double active_ree_y2_sum = 0.0;
    double active_ree_z2_sum = 0.0;
    double active_stretch_sum = 0.0;
    double active_p2z_sum = 0.0;

    long long contour_segments = 0;
    double contour_length = 0.0;
    double active_contour_length = 0.0;
    double defect_contour_length = 0.0;
    double segment_p2x_sum = 0.0;
    double segment_p2y_sum = 0.0;
    double segment_p2z_sum = 0.0;
    long long active_contour_segments = 0;
    double active_segment_p2x_sum = 0.0;
    double active_segment_p2y_sum = 0.0;
    double active_segment_p2z_sum = 0.0;

    long long z1_kinks = 0;
    long long z1_segments = 0;
    double z1_primitive_length = 0.0;
    double z1_segment_p2z_sum = 0.0;
};

struct Z1Point {
    Vec3 position;
    double contour_index = 0.0;
    int kink = 0;
};

struct Z1Result {
    Box box;
    std::vector<std::vector<Z1Point>> chains;
};

struct DumpFrame {
    long long timestep = 0;
    Box box;
    std::vector<Vec3> unwrapped;
};

struct LinearFit {
    long long points = 0;
    double slope = std::numeric_limits<double>::quiet_NaN();
    double intercept = std::numeric_limits<double>::quiet_NaN();
    double r_squared = std::numeric_limits<double>::quiet_NaN();
};

LinearFit linear_fit(const std::vector<double> &time,
                     const std::vector<double> &value,
                     double start_time,
                     double end_time =
                         std::numeric_limits<double>::infinity()) {
    LinearFit result;
    double sum_x = 0.0, sum_y = 0.0;
    for (std::size_t index = 0; index < time.size(); ++index) {
        if (time[index] < start_time || time[index] > end_time ||
            !std::isfinite(value[index])) continue;
        ++result.points;
        sum_x += time[index];
        sum_y += value[index];
    }
    if (result.points < 2) return result;
    const double mean_x = sum_x / result.points;
    const double mean_y = sum_y / result.points;
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t index = 0; index < time.size(); ++index) {
        if (time[index] < start_time || time[index] > end_time ||
            !std::isfinite(value[index])) continue;
        const double dx = time[index] - mean_x;
        const double dy = value[index] - mean_y;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (!(sxx > 0.0)) return result;
    result.slope = sxy / sxx;
    result.intercept = mean_y - result.slope * mean_x;
    double residual = 0.0;
    for (std::size_t index = 0; index < time.size(); ++index) {
        if (time[index] < start_time || time[index] > end_time ||
            !std::isfinite(value[index])) continue;
        const double difference = value[index] -
            (result.intercept + result.slope * time[index]);
        residual += difference * difference;
    }
    if (syy > 0.0) result.r_squared = 1.0 - residual / syy;
    return result;
}

std::vector<std::string> words(const std::string &line) {
    std::vector<std::string> result;
    std::istringstream fields(line);
    std::string value;
    while (fields >> value) result.push_back(value);
    return result;
}

class DumpReader {
  public:
    explicit DumpReader(const std::string &path) : input_(path), path_(path) {
        if (!input_) throw std::runtime_error("cannot open trajectory: " + path);
    }

    bool next(DumpFrame &frame, long long expected_atoms) {
        std::string line;
        while (std::getline(input_, line) && trim(line).empty()) {}
        if (!input_) return false;
        if (trim(line) != "ITEM: TIMESTEP")
            throw std::runtime_error("expected ITEM: TIMESTEP in " + path_);
        if (!std::getline(input_, line)) throw std::runtime_error("missing timestep");
        frame.timestep = std::stoll(trim(line));
        require("ITEM: NUMBER OF ATOMS");
        if (!std::getline(input_, line)) throw std::runtime_error("missing atom count");
        const long long atom_count = std::stoll(trim(line));
        if (atom_count != expected_atoms)
            throw std::runtime_error("trajectory atom count differs from topology data");
        if (!std::getline(input_, line) || !begins_with(trim(line), "ITEM: BOX BOUNDS"))
            throw std::runtime_error("missing BOX BOUNDS");
        frame.box = {};
        read_bounds(frame.box.xlo, frame.box.xhi);
        read_bounds(frame.box.ylo, frame.box.yhi);
        read_bounds(frame.box.zlo, frame.box.zhi);
        frame.box.have_x = frame.box.have_y = frame.box.have_z = true;
        if (!std::getline(input_, line) || !begins_with(trim(line), "ITEM: ATOMS"))
            throw std::runtime_error("missing ATOMS header");
        std::vector<std::string> columns = words(line);
        columns.erase(columns.begin(), columns.begin() + 2);
        const auto column = [&](const std::string &name) {
            const auto found = std::find(columns.begin(), columns.end(), name);
            return found == columns.end() ? -1 :
                static_cast<int>(found - columns.begin());
        };
        const int id_column = column("id");
        const int xu_column = column("xu"), yu_column = column("yu"),
                  zu_column = column("zu");
        const int x_column = column("x"), y_column = column("y"),
                  z_column = column("z");
        const int ix_column = column("ix"), iy_column = column("iy"),
                  iz_column = column("iz");
        const bool have_unwrapped =
            xu_column >= 0 && yu_column >= 0 && zu_column >= 0;
        const bool have_wrapped_images = x_column >= 0 && y_column >= 0 &&
            z_column >= 0 && ix_column >= 0 && iy_column >= 0 && iz_column >= 0;
        if (id_column < 0 || (!have_unwrapped && !have_wrapped_images))
            throw std::runtime_error(
                "trajectory needs id plus xu/yu/zu or x/y/z/ix/iy/iz columns");
        frame.unwrapped.assign(static_cast<std::size_t>(atom_count + 1), {
            std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0});
        for (long long row = 0; row < atom_count; ++row) {
            if (!std::getline(input_, line)) throw std::runtime_error("truncated ATOMS block");
            const std::vector<std::string> values = words(line);
            if (values.size() < columns.size()) throw std::runtime_error("short atom row");
            const long long id = std::stoll(values[static_cast<std::size_t>(id_column)]);
            if (id < 1 || id > atom_count ||
                std::isfinite(frame.unwrapped[static_cast<std::size_t>(id)].x))
                throw std::runtime_error("invalid or duplicate trajectory atom ID");
            Vec3 position;
            if (have_unwrapped) {
                position = {
                    std::stod(values[static_cast<std::size_t>(xu_column)]),
                    std::stod(values[static_cast<std::size_t>(yu_column)]),
                    std::stod(values[static_cast<std::size_t>(zu_column)])};
            } else {
                position = {
                    std::stod(values[static_cast<std::size_t>(x_column)]) +
                        std::stoll(values[static_cast<std::size_t>(ix_column)]) * frame.box.lx(),
                    std::stod(values[static_cast<std::size_t>(y_column)]) +
                        std::stoll(values[static_cast<std::size_t>(iy_column)]) * frame.box.ly(),
                    std::stod(values[static_cast<std::size_t>(z_column)]) +
                        std::stoll(values[static_cast<std::size_t>(iz_column)]) * frame.box.lz()};
            }
            frame.unwrapped[static_cast<std::size_t>(id)] = position;
        }
        return true;
    }

  private:
    void require(const std::string &expected) {
        std::string line;
        if (!std::getline(input_, line) || trim(line) != expected)
            throw std::runtime_error("expected " + expected + " in " + path_);
    }
    void read_bounds(double &lo, double &hi) {
        std::string line;
        if (!std::getline(input_, line)) throw std::runtime_error("truncated BOX BOUNDS");
        std::istringstream fields(line);
        double tilt = 0.0;
        if (!(fields >> lo >> hi)) throw std::runtime_error("invalid box bounds");
        if (fields >> tilt)
            throw std::runtime_error("triclinic trajectories are not currently supported");
    }
    std::ifstream input_;
    std::string path_;
};

int bin_index(double z, const Box &box, int bins, bool periodic_z) {
    if (periodic_z) z = wrap_position({0.0, 0.0, z}, box, true).z;
    int index = static_cast<int>(std::floor((z - box.zlo) / box.lz() * bins));
    if (index < 0) index = 0;
    if (index >= bins) index = bins - 1;
    return index;
}

Vec3 path_midpoint(
    const EffectiveStrand &strand, const DataFile &data,
    const ModelInfo &info) {
    const auto positions = unwrapped_path(strand.atoms, data, info);
    const double half = 0.5 * strand.contour_length;
    double traversed = 0.0;
    for (std::size_t i = 1; i < positions.size(); ++i) {
        const Vec3 segment = positions[i] - positions[i - 1];
        const double length = norm(segment);
        if (traversed + length >= half && length > 0.0) {
            const Vec3 point = positions[i - 1] +
                ((half - traversed) / length) * segment;
            return wrap_position(point, data.box, info.periodic_z());
        }
        traversed += length;
    }
    return wrap_position(positions.back(), data.box, info.periodic_z());
}

Vec3 molecule_center(
    long long molecule, const DataFile &data, const ModelInfo &info) {
    const auto &atoms = data.molecule_atoms[static_cast<std::size_t>(molecule)];
    const Vec3 reference = data.atoms[static_cast<std::size_t>(atoms.front())].position;
    Vec3 sum = reference;
    for (std::size_t i = 1; i < atoms.size(); ++i) {
        const Vec3 position = data.atoms[static_cast<std::size_t>(atoms[i])].position;
        sum += reference + minimum_image_vector(
            position - reference, data.box, info.periodic_z());
    }
    return wrap_position((1.0 / atoms.size()) * sum, data.box, info.periodic_z());
}

double mean_or_nan(double sum, long long count) {
    return count > 0 ? sum / static_cast<double>(count) :
        std::numeric_limits<double>::quiet_NaN();
}

double ratio_or_nan(double numerator, double denominator) {
    return denominator != 0.0 ? numerator / denominator :
        std::numeric_limits<double>::quiet_NaN();
}

void add_profile(BinProfile &target, const BinProfile &source) {
    target.topology.crosslink_bonds += source.topology.crosslink_bonds;
    target.topology.strand_crosslink_bonds +=
        source.topology.strand_crosslink_bonds;
    target.topology.moderator_crosslink_bonds +=
        source.topology.moderator_crosslink_bonds;
    target.topology.junctions += source.topology.junctions;
    target.topology.active_strands += source.topology.active_strands;
    target.topology.dangling_ends += source.topology.dangling_ends;
    target.topology.dangling_loops += source.topology.dangling_loops;
    target.topology.self_loops += source.topology.self_loops;
    target.topology.isolated_parents += source.topology.isolated_parents;
    for (int component = 0; component < kComponentCount; ++component) {
        target.beads[static_cast<std::size_t>(component)] +=
            source.beads[static_cast<std::size_t>(component)];
        target.mass[static_cast<std::size_t>(component)] +=
            source.mass[static_cast<std::size_t>(component)];
    }
    for (std::size_t component = 0; component < 3; ++component) {
        target.functional_total[component] += source.functional_total[component];
        target.functional_reacted[component] += source.functional_reacted[component];
    }
#define ADD_FIELD(field) target.field += source.field
    ADD_FIELD(strand_markers);
    ADD_FIELD(strand_orientation_count);
    ADD_FIELD(strand_ree_x2_sum); ADD_FIELD(strand_ree_y2_sum);
    ADD_FIELD(strand_ree_z2_sum); ADD_FIELD(strand_stretch_sum);
    ADD_FIELD(strand_p2z_sum);
    ADD_FIELD(active_orientation_count);
    ADD_FIELD(active_ree_x2_sum); ADD_FIELD(active_ree_y2_sum);
    ADD_FIELD(active_ree_z2_sum); ADD_FIELD(active_stretch_sum);
    ADD_FIELD(active_p2z_sum);
    ADD_FIELD(contour_segments); ADD_FIELD(contour_length);
    ADD_FIELD(active_contour_length); ADD_FIELD(defect_contour_length);
    ADD_FIELD(segment_p2x_sum); ADD_FIELD(segment_p2y_sum);
    ADD_FIELD(segment_p2z_sum); ADD_FIELD(active_contour_segments);
    ADD_FIELD(active_segment_p2x_sum); ADD_FIELD(active_segment_p2y_sum);
    ADD_FIELD(active_segment_p2z_sum);
    ADD_FIELD(z1_kinks); ADD_FIELD(z1_segments);
    ADD_FIELD(z1_primitive_length); ADD_FIELD(z1_segment_p2z_sum);
#undef ADD_FIELD
}

Z1Result read_z1_sp(const std::string &path, const DataFile &data) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Z1+ result: " + path);
    long long chain_count = 0;
    Z1Result result;
    if (!(input >> chain_count) || chain_count < 0 ||
        !(input >> result.box.xhi >> result.box.yhi >> result.box.zhi))
        throw std::runtime_error("invalid Z1+SP header in " + path);
    result.box.have_x = result.box.have_y = result.box.have_z = true;
    const auto close_length = [](double first, double second) {
        return std::fabs(first - second) <=
            1.0e-6 * std::max({1.0, std::fabs(first), std::fabs(second)});
    };
    if (!close_length(result.box.lx(), data.box.lx()) ||
        !close_length(result.box.ly(), data.box.ly()) ||
        !close_length(result.box.lz(), data.box.lz()))
        throw std::runtime_error(
            "Z1+SP box lengths differ from the LAMMPS snapshot; "
            "use unscaled Z1+ coordinates for a physical z profile");
    result.chains.reserve(static_cast<std::size_t>(chain_count));
    for (long long chain = 0; chain < chain_count; ++chain) {
        long long point_count = 0;
        if (!(input >> point_count) || point_count < 2)
            throw std::runtime_error("invalid Z1+ primitive path point count");
        std::vector<Z1Point> points(static_cast<std::size_t>(point_count));
        for (Z1Point &point : points)
            if (!(input >> point.position.x >> point.position.y >> point.position.z >>
                  point.contour_index >> point.kink))
                throw std::runtime_error("truncated Z1+ primitive path data");
        result.chains.push_back(std::move(points));
    }
    std::string trailing;
    if (input >> trailing)
        throw std::runtime_error("unexpected trailing data in Z1+SP result");
    return result;
}

std::set<long long> reacted_functional_atoms(
    const DataFile &data, const ModelInfo &info) {
    std::set<long long> reacted;
    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        reacted.insert(bond.first);
        reacted.insert(bond.second);
    }
    return reacted;
}

std::vector<BinProfile> static_profile(
    const DataFile &data, const ModelInfo &info,
    const ReducedNetwork &network, int bins,
    const Z1Result *z1_result) {
    std::vector<BinProfile> profile(static_cast<std::size_t>(bins));
    const std::set<long long> reacted = reacted_functional_atoms(data, info);

    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        const int component = component_for_molecule(atom.molecule, info);
        BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
            atom.position.z, data.box, bins, info.periodic_z()))];
        ++bin.beads[static_cast<std::size_t>(component)];
        bin.mass[static_cast<std::size_t>(component)] += atom.mass;
    }

    const auto add_declared_sites = [&](int component,
                                        const std::vector<long long> &sites,
                                        std::size_t profile_component) {
        const ComponentInfo &entry = info.components[static_cast<std::size_t>(component)];
        const std::set<long long> declared(sites.begin(), sites.end());
        for (long long molecule = entry.molecule_start;
             molecule > 0 && molecule <= entry.molecule_end; ++molecule) {
            const auto &atoms = data.molecule_atoms[static_cast<std::size_t>(molecule)];
            for (long long atom_id : atoms) {
                if (!declared.count(local_rank(atom_id, atoms))) continue;
                BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
                    data.atoms[static_cast<std::size_t>(atom_id)].position.z,
                    data.box, bins, info.periodic_z()))];
                ++bin.functional_total[profile_component];
                if (reacted.count(atom_id))
                    ++bin.functional_reacted[profile_component];
            }
        }
    };
    add_declared_sites(kStrand, info.reactive_bead_sites, 0);
    add_declared_sites(kCrosslinker, info.crosslinker_reactive_bead_sites, 1);
    const ComponentInfo &moderators = info.components[kModerator];
    for (long long molecule = moderators.molecule_start;
         molecule > 0 && molecule <= moderators.molecule_end; ++molecule) {
        for (long long atom_id :
             data.molecule_atoms[static_cast<std::size_t>(molecule)]) {
            const Atom &atom = data.atoms[static_cast<std::size_t>(atom_id)];
            // The implemented five-bead moderator has a neutral center and
            // four type-2 functional arms.
            if (atom.type != 2) continue;
            BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
                atom.position.z, data.box, bins, info.periodic_z()))];
            ++bin.functional_total[2];
            if (reacted.count(atom_id)) ++bin.functional_reacted[2];
        }
    }

    for (const Bond &bond : data.bonds) {
        if (bond.type != info.crosslink_bond_type) continue;
        const Atom &first = data.atoms[static_cast<std::size_t>(bond.first)];
        const Atom &second = data.atoms[static_cast<std::size_t>(bond.second)];
        const Vec3 delta = minimum_image_vector(
            second.position - first.position, data.box, info.periodic_z());
        const Vec3 midpoint = wrap_position(
            first.position + 0.5 * delta, data.box, info.periodic_z());
        BinCounts &bin = profile[static_cast<std::size_t>(
            bin_index(midpoint.z, data.box, bins, info.periodic_z()))].topology;
        ++bin.crosslink_bonds;
        const int first_component = component_for_molecule(first.molecule, info);
        const int second_component = component_for_molecule(second.molecule, info);
        if (first_component == kStrand || second_component == kStrand)
            ++bin.strand_crosslink_bonds;
        if (first_component == kModerator || second_component == kModerator)
            ++bin.moderator_crosslink_bonds;
    }
    for (std::size_t node = 1; node < network.nodes.size(); ++node)
        ++profile[static_cast<std::size_t>(bin_index(
            network.nodes[node].position.z, data.box, bins,
            info.periodic_z()))].topology.junctions;
    for (const EffectiveStrand &strand : network.strands) {
        Vec3 marker = path_midpoint(strand, data, info);
        BinProfile &marker_bin = profile[static_cast<std::size_t>(bin_index(
            marker.z, data.box, bins, info.periodic_z()))];
        BinCounts &bin = marker_bin.topology;
        ++marker_bin.strand_markers;
        if (strand.end_to_end_length > 0.0) {
            const double inverse_ree2 =
                1.0 / (strand.end_to_end_length * strand.end_to_end_length);
            marker_bin.strand_ree_x2_sum += strand.end_to_end.x * strand.end_to_end.x;
            marker_bin.strand_ree_y2_sum += strand.end_to_end.y * strand.end_to_end.y;
            marker_bin.strand_ree_z2_sum += strand.end_to_end.z * strand.end_to_end.z;
            marker_bin.strand_p2z_sum +=
                0.5 * (3.0 * strand.end_to_end.z * strand.end_to_end.z * inverse_ree2 - 1.0);
            marker_bin.strand_stretch_sum +=
                ratio_or_nan(strand.end_to_end_length, strand.contour_length);
            ++marker_bin.strand_orientation_count;
        }
        if (strand.status == "active") {
            ++bin.active_strands;
            if (strand.end_to_end_length > 0.0) {
                const double inverse_ree2 =
                    1.0 / (strand.end_to_end_length * strand.end_to_end_length);
                marker_bin.active_ree_x2_sum += strand.end_to_end.x * strand.end_to_end.x;
                marker_bin.active_ree_y2_sum += strand.end_to_end.y * strand.end_to_end.y;
                marker_bin.active_ree_z2_sum += strand.end_to_end.z * strand.end_to_end.z;
                marker_bin.active_p2z_sum +=
                    0.5 * (3.0 * strand.end_to_end.z * strand.end_to_end.z * inverse_ree2 - 1.0);
                marker_bin.active_stretch_sum +=
                    ratio_or_nan(strand.end_to_end_length, strand.contour_length);
                ++marker_bin.active_orientation_count;
            }
        } else if (strand.status == "dangling_loop") ++bin.dangling_loops;
        else if (strand.status == "self_loop") ++bin.self_loops;
        else if (strand.status == "dangling") {
            const long long free_atom = strand.first_node == 0
                ? strand.first_atom : strand.second_atom;
            const double z = data.atoms[static_cast<std::size_t>(free_atom)].position.z;
            ++profile[static_cast<std::size_t>(bin_index(
                z, data.box, bins, info.periodic_z()))].topology.dangling_ends;
        }

        const std::vector<Vec3> positions = unwrapped_path(strand.atoms, data, info);
        for (std::size_t point = 1; point < positions.size(); ++point) {
            const Vec3 segment = positions[point] - positions[point - 1];
            const double length = norm(segment);
            if (!(length > 0.0)) continue;
            const Vec3 midpoint = wrap_position(
                positions[point - 1] + 0.5 * segment,
                data.box, info.periodic_z());
            BinProfile &segment_bin = profile[static_cast<std::size_t>(bin_index(
                midpoint.z, data.box, bins, info.periodic_z()))];
            const double inverse_length2 = 1.0 / (length * length);
            const double p2x = 0.5 * (3.0 * segment.x * segment.x * inverse_length2 - 1.0);
            const double p2y = 0.5 * (3.0 * segment.y * segment.y * inverse_length2 - 1.0);
            const double p2z = 0.5 * (3.0 * segment.z * segment.z * inverse_length2 - 1.0);
            ++segment_bin.contour_segments;
            segment_bin.contour_length += length;
            segment_bin.segment_p2x_sum += p2x;
            segment_bin.segment_p2y_sum += p2y;
            segment_bin.segment_p2z_sum += p2z;
            if (strand.status == "active") {
                segment_bin.active_contour_length += length;
                ++segment_bin.active_contour_segments;
                segment_bin.active_segment_p2x_sum += p2x;
                segment_bin.active_segment_p2y_sum += p2y;
                segment_bin.active_segment_p2z_sum += p2z;
            } else {
                segment_bin.defect_contour_length += length;
            }
        }
    }
    for (const ParentRecord &parent : network.parents) {
        if (parent.state != "isolated" && parent.state != "isolated_ring" &&
            parent.state != "isolated_star" && parent.state != "isolated_grafted") continue;
        const Vec3 center = molecule_center(parent.molecule, data, info);
        ++profile[static_cast<std::size_t>(bin_index(
            center.z, data.box, bins, info.periodic_z()))].topology.isolated_parents;
    }

    if (z1_result != nullptr) {
        for (const auto &chain : z1_result->chains) {
            for (const Z1Point &point : chain) {
                if (point.kink == 0) continue;
                ++profile[static_cast<std::size_t>(bin_index(
                    point.position.z, data.box, bins,
                    info.periodic_z()))].z1_kinks;
            }
            for (std::size_t point = 1; point < chain.size(); ++point) {
                const Vec3 segment =
                    chain[point].position - chain[point - 1].position;
                const double length = norm(segment);
                if (!(length > 0.0)) continue;
                Vec3 midpoint = chain[point - 1].position + 0.5 * segment;
                midpoint = wrap_position(midpoint, data.box, info.periodic_z());
                BinProfile &bin = profile[static_cast<std::size_t>(bin_index(
                    midpoint.z, data.box, bins, info.periodic_z()))];
                ++bin.z1_segments;
                bin.z1_primitive_length += length;
                bin.z1_segment_p2z_sum += 0.5 *
                    (3.0 * segment.z * segment.z / (length * length) - 1.0);
            }
        }
    }
    return profile;
}

struct SlidingProfile {
    SlidingGrid grid;
    std::vector<int> center_indices;
    std::vector<BinProfile> windows;
};

SlidingProfile make_sliding_profile(
    const std::vector<BinProfile> &fine_profile, const SlidingGrid &grid,
    bool periodic_z) {
    SlidingProfile result;
    result.grid = grid;
    const int first_center = periodic_z ? 0 : grid.half_window_bins;
    const int final_center = periodic_z
        ? grid.fine_bins : grid.fine_bins - grid.half_window_bins;
    result.center_indices.reserve(
        static_cast<std::size_t>(std::max(0, final_center - first_center)));
    result.windows.reserve(result.center_indices.capacity());
    for (int center = first_center; center < final_center; ++center) {
        BinProfile window;
        for (int offset = -grid.half_window_bins;
             offset <= grid.half_window_bins; ++offset) {
            int fine_bin = center + offset;
            if (periodic_z) {
                fine_bin %= grid.fine_bins;
                if (fine_bin < 0) fine_bin += grid.fine_bins;
            }
            add_profile(
                window, fine_profile[static_cast<std::size_t>(fine_bin)]);
        }
        result.center_indices.push_back(center);
        result.windows.push_back(std::move(window));
    }
    return result;
}

void write_profile(
    const std::filesystem::path &path, const std::vector<BinProfile> &profile,
    const DataFile &data, const ModelInfo &info, bool z1_available,
    const SlidingProfile *sliding = nullptr) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    if (sliding == nullptr)
        out << "bin\tzlo_A\tzhi_A\tz_center_A\tvolume_A3";
    else
        out << "window\tz_center_A\tzlo_A\tzhi_A\tstep_A"
            << "\twindow_width_A\tvolume_A3";
    out << "\tcrosslink_bonds"
        << "\tstrand_crosslink_bonds\tmoderator_crosslink_bonds\tjunctions"
        << "\tactive_strands\tdangling_ends\tdangling_loops\tself_loops"
        << "\tisolated_parents\tcrosslink_density_A-3\tdefect_density_A-3"
        << "\tmaterial_z_normalized\tdistance_from_nearest_box_wall_A"
        << "\tdistance_from_effective_wall_A"
        << "\tbeads_total\tbeads_strand\tbeads_crosslinker\tbeads_moderator\tbeads_filler"
        << "\tmass_density_total_g_cm-3\tmass_density_strand_g_cm-3"
        << "\tmass_density_crosslinker_g_cm-3\tmass_density_moderator_g_cm-3"
        << "\tmass_density_filler_g_cm-3"
        << "\tfunctional_sites_total\tfunctional_sites_reacted\tlocal_conversion"
        << "\tstrand_sites_total\tstrand_sites_reacted\tstrand_local_conversion"
        << "\tcrosslinker_sites_total\tcrosslinker_sites_reacted"
        << "\tcrosslinker_local_conversion\tmoderator_sites_total"
        << "\tmoderator_sites_reacted\tmoderator_local_conversion"
        << "\tjunction_density_A-3\tactive_strand_density_A-3"
        << "\tdangling_end_density_A-3\tdangling_loop_density_A-3"
        << "\tself_loop_density_A-3\tisolated_parent_density_A-3"
        << "\tcontour_length_A\tactive_contour_length_A\tdefect_contour_length_A"
        << "\tcontour_length_density_A-2\tactive_contour_length_density_A-2"
        << "\tdefect_contour_length_density_A-2\tcontour_segments"
        << "\tsegment_P2_x\tsegment_P2_y\tsegment_P2_z"
        << "\tactive_contour_segments\tactive_segment_P2_x"
        << "\tactive_segment_P2_y\tactive_segment_P2_z"
        << "\tstrand_markers\tstrand_Ree_x2_mean_A2\tstrand_Ree_y2_mean_A2"
        << "\tstrand_Ree_z2_mean_A2\tstrand_Ree_parallel2_mean_A2"
        << "\tstrand_Ree_over_Lc_mean\tstrand_P2_z"
        << "\tactive_Ree_x2_mean_A2\tactive_Ree_y2_mean_A2"
        << "\tactive_Ree_z2_mean_A2\tactive_Ree_parallel2_mean_A2"
        << "\tactive_Ree_over_Lc_mean\tactive_strand_P2_z"
        << "\tz1_data_available\tz1_kinks\tz1_kink_density_A-3"
        << "\tz1_primitive_segments\tz1_primitive_length_A"
        << "\tz1_primitive_length_density_A-2\tz1_primitive_P2_z\n";
    const double width = sliding == nullptr
        ? data.box.lz() / profile.size()
        : sliding->grid.window_width_A;
    const double volume = data.box.lx() * data.box.ly() * width;
    const double wall_cutoff = info.geometry == "film"
        ? info.film_wall_cutoff_per_side_angstrom : 0.0;
    double material_thickness = info.geometry == "film"
        ? info.film_thickness_angstrom : data.box.lz();
    if (!(material_thickness > 0.0))
        material_thickness = data.box.lz() - 2.0 * wall_cutoff;
    out << std::setprecision(12);
    for (std::size_t bin = 0; bin < profile.size(); ++bin) {
        const BinProfile &entry = profile[bin];
        const BinCounts &counts = entry.topology;
        const double center = sliding == nullptr
            ? data.box.zlo + (bin + 0.5) * width
            : data.box.zlo +
                (sliding->center_indices[bin] + 0.5) * sliding->grid.step_A;
        const double zlo = center - 0.5 * width;
        const double distance_box = std::min(
            center - data.box.zlo, data.box.zhi - center);
        const double material_z =
            (center - data.box.zlo - wall_cutoff) / material_thickness;
        const long long defects = counts.dangling_ends + counts.dangling_loops +
            counts.self_loops + counts.isolated_parents;
        long long beads_total = 0;
        double mass_total = 0.0;
        long long sites_total = 0, sites_reacted = 0;
        for (int component = 0; component < kComponentCount; ++component) {
            beads_total += entry.beads[static_cast<std::size_t>(component)];
            mass_total += entry.mass[static_cast<std::size_t>(component)];
        }
        for (std::size_t component = 0; component < 3; ++component) {
            sites_total += entry.functional_total[component];
            sites_reacted += entry.functional_reacted[component];
        }
        const auto density = [&](double value) { return value / volume; };
        const auto mass_density = [&](double mass) {
            return mass / (kAvogadroAngstrom * volume);
        };
        const double z1_kink_density = z1_available
            ? density(entry.z1_kinks) : std::numeric_limits<double>::quiet_NaN();
        const double z1_length_density = z1_available
            ? density(entry.z1_primitive_length) :
              std::numeric_limits<double>::quiet_NaN();
        const double z1_p2z = z1_available
            ? mean_or_nan(entry.z1_segment_p2z_sum, entry.z1_segments) :
              std::numeric_limits<double>::quiet_NaN();
        if (sliding == nullptr)
            out << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
                << center << '\t' << volume;
        else
            out << bin + 1 << '\t' << center << '\t' << zlo << '\t'
                << zlo + width << '\t' << sliding->grid.step_A << '\t'
                << width << '\t' << volume;
        out << '\t' << counts.crosslink_bonds << '\t'
            << counts.strand_crosslink_bonds
            << '\t' << counts.moderator_crosslink_bonds << '\t' << counts.junctions
            << '\t' << counts.active_strands << '\t' << counts.dangling_ends
            << '\t' << counts.dangling_loops << '\t' << counts.self_loops
            << '\t' << counts.isolated_parents << '\t'
            << density(counts.crosslink_bonds) << '\t' << density(defects) << '\t'
            << material_z << '\t' << distance_box << '\t'
            << distance_box - wall_cutoff << '\t'
            << beads_total;
        for (long long value : entry.beads) out << '\t' << value;
        out << '\t' << mass_density(mass_total);
        for (double mass : entry.mass) out << '\t' << mass_density(mass);
        out << '\t' << sites_total << '\t' << sites_reacted << '\t'
            << ratio_or_nan(sites_reacted, sites_total);
        for (std::size_t component = 0; component < 3; ++component)
            out << '\t' << entry.functional_total[component] << '\t'
                << entry.functional_reacted[component] << '\t'
                << ratio_or_nan(entry.functional_reacted[component],
                                entry.functional_total[component]);
        out << '\t' << density(counts.junctions)
            << '\t' << density(counts.active_strands)
            << '\t' << density(counts.dangling_ends)
            << '\t' << density(counts.dangling_loops)
            << '\t' << density(counts.self_loops)
            << '\t' << density(counts.isolated_parents)
            << '\t' << entry.contour_length
            << '\t' << entry.active_contour_length
            << '\t' << entry.defect_contour_length
            << '\t' << density(entry.contour_length)
            << '\t' << density(entry.active_contour_length)
            << '\t' << density(entry.defect_contour_length)
            << '\t' << entry.contour_segments
            << '\t' << mean_or_nan(entry.segment_p2x_sum, entry.contour_segments)
            << '\t' << mean_or_nan(entry.segment_p2y_sum, entry.contour_segments)
            << '\t' << mean_or_nan(entry.segment_p2z_sum, entry.contour_segments)
            << '\t' << entry.active_contour_segments
            << '\t' << mean_or_nan(entry.active_segment_p2x_sum,
                                    entry.active_contour_segments)
            << '\t' << mean_or_nan(entry.active_segment_p2y_sum,
                                    entry.active_contour_segments)
            << '\t' << mean_or_nan(entry.active_segment_p2z_sum,
                                    entry.active_contour_segments)
            << '\t' << entry.strand_markers
            << '\t' << mean_or_nan(entry.strand_ree_x2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_ree_y2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_ree_z2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_ree_x2_sum + entry.strand_ree_y2_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_stretch_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.strand_p2z_sum,
                                    entry.strand_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_x2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_y2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_z2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_ree_x2_sum + entry.active_ree_y2_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_stretch_sum,
                                    entry.active_orientation_count)
            << '\t' << mean_or_nan(entry.active_p2z_sum,
                                    entry.active_orientation_count)
            << '\t' << (z1_available ? 1 : 0)
            << '\t' << entry.z1_kinks << '\t' << z1_kink_density
            << '\t' << entry.z1_segments << '\t' << entry.z1_primitive_length
            << '\t' << z1_length_density << '\t' << z1_p2z << '\n';
    }
}

void write_folded_profile(
    const std::filesystem::path &path, const std::vector<BinProfile> &profile,
    const DataFile &data, const ModelInfo &info, bool z1_available,
    const SlidingProfile *sliding = nullptr) {
    const std::size_t folded_bins = (profile.size() + 1) / 2;
    std::vector<BinProfile> folded(folded_bins);
    std::vector<int> multiplicity(folded_bins, 0);
    for (std::size_t bin = 0; bin < profile.size(); ++bin) {
        const std::size_t folded_bin = std::min(bin, profile.size() - 1 - bin);
        add_profile(folded[folded_bin], profile[bin]);
        ++multiplicity[folded_bin];
    }
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "folded_bin\tdistance_box_lo_A\tdistance_box_hi_A"
        << "\tdistance_box_center_A";
    if (sliding != nullptr)
        out << "\tstep_A\twindow_width_A";
    out << "\tdistance_effective_wall_center_A"
        << "\tnormalized_distance_effective_wall\tpaired_full_bins\tvolume_A3"
        << "\tmass_density_total_g_cm-3\tcrosslink_density_A-3"
        << "\tjunction_density_A-3\tactive_strand_density_A-3"
        << "\tdefect_density_A-3\tlocal_conversion"
        << "\tactive_contour_length_density_A-2\tdefect_contour_length_density_A-2"
        << "\tsegment_P2_z\tactive_segment_P2_z\tactive_strand_P2_z"
        << "\tactive_Ree_parallel2_mean_A2\tactive_Ree_z2_mean_A2"
        << "\tz1_data_available\tz1_kinks\tz1_kink_density_A-3"
        << "\tz1_primitive_length_density_A-2\tz1_primitive_P2_z\n";
    const double width = sliding == nullptr
        ? data.box.lz() / profile.size()
        : sliding->grid.window_width_A;
    const double area = data.box.lx() * data.box.ly();
    const double wall_cutoff = info.geometry == "film"
        ? info.film_wall_cutoff_per_side_angstrom : 0.0;
    double material_thickness = info.geometry == "film"
        ? info.film_thickness_angstrom : data.box.lz();
    if (!(material_thickness > 0.0))
        material_thickness = data.box.lz() - 2.0 * wall_cutoff;
    out << std::setprecision(12);
    for (std::size_t bin = 0; bin < folded.size(); ++bin) {
        const BinProfile &entry = folded[bin];
        const double center = sliding == nullptr
            ? (bin + 0.5) * width
            : std::min(
                data.box.zlo +
                    (sliding->center_indices[bin] + 0.5) *
                        sliding->grid.step_A - data.box.zlo,
                data.box.zhi -
                    (data.box.zlo +
                     (sliding->center_indices[bin] + 0.5) *
                         sliding->grid.step_A));
        const double distance_lo = std::max(0.0, center - 0.5 * width);
        const double distance_hi = std::min(
            center + 0.5 * width, 0.5 * data.box.lz());
        const double volume = area * width * multiplicity[bin];
        double mass = 0.0;
        long long total_sites = 0, reacted_sites = 0;
        for (double value : entry.mass) mass += value;
        for (std::size_t component = 0; component < 3; ++component) {
            total_sites += entry.functional_total[component];
            reacted_sites += entry.functional_reacted[component];
        }
        const long long defects = entry.topology.dangling_ends +
            entry.topology.dangling_loops + entry.topology.self_loops +
            entry.topology.isolated_parents;
        const auto density = [&](double value) { return value / volume; };
        out << bin + 1 << '\t' << distance_lo << '\t' << distance_hi << '\t'
            << center;
        if (sliding != nullptr)
            out << '\t' << sliding->grid.step_A << '\t' << width;
        out << '\t' << center - wall_cutoff << '\t'
            << 2.0 * (center - wall_cutoff) / material_thickness << '\t'
            << multiplicity[bin] << '\t' << volume << '\t'
            << mass / (kAvogadroAngstrom * volume) << '\t'
            << density(entry.topology.crosslink_bonds) << '\t'
            << density(entry.topology.junctions) << '\t'
            << density(entry.topology.active_strands) << '\t'
            << density(defects) << '\t'
            << ratio_or_nan(reacted_sites, total_sites) << '\t'
            << density(entry.active_contour_length) << '\t'
            << density(entry.defect_contour_length) << '\t'
            << mean_or_nan(entry.segment_p2z_sum, entry.contour_segments) << '\t'
            << mean_or_nan(entry.active_segment_p2z_sum,
                            entry.active_contour_segments) << '\t'
            << mean_or_nan(entry.active_p2z_sum,
                            entry.active_orientation_count) << '\t'
            << mean_or_nan(entry.active_ree_x2_sum + entry.active_ree_y2_sum,
                            entry.active_orientation_count) << '\t'
            << mean_or_nan(entry.active_ree_z2_sum,
                            entry.active_orientation_count) << '\t'
            << (z1_available ? 1 : 0) << '\t' << entry.z1_kinks << '\t'
            << (z1_available ? density(entry.z1_kinks) :
                std::numeric_limits<double>::quiet_NaN()) << '\t'
            << (z1_available ? density(entry.z1_primitive_length) :
                std::numeric_limits<double>::quiet_NaN()) << '\t'
            << (z1_available ? mean_or_nan(entry.z1_segment_p2z_sum,
                                           entry.z1_segments) :
                std::numeric_limits<double>::quiet_NaN()) << '\n';
    }
}

struct RegionProfile {
    BinProfile profile;
    double volume = 0.0;
    long long bins = 0;
};

void write_profile_summary(
    const std::filesystem::path &path, const std::vector<BinProfile> &profile,
    const DataFile &data, const ModelInfo &info, const Options &options,
    bool z1_available) {
    const double width = data.box.lz() / profile.size();
    const double bin_volume = data.box.lx() * data.box.ly() * width;
    const double wall_cutoff = info.geometry == "film"
        ? info.film_wall_cutoff_per_side_angstrom : 0.0;
    double material_thickness = info.geometry == "film"
        ? info.film_thickness_angstrom : data.box.lz();
    if (!(material_thickness > 0.0))
        material_thickness = data.box.lz() - 2.0 * wall_cutoff;
    RegionProfile lower, upper, core;
    const auto include = [&](RegionProfile &region, const BinProfile &entry) {
        add_profile(region.profile, entry);
        region.volume += bin_volume;
        ++region.bins;
    };
    for (std::size_t bin = 0; bin < profile.size(); ++bin) {
        const double center = data.box.zlo + (bin + 0.5) * width;
        const double normalized =
            (center - data.box.zlo - wall_cutoff) / material_thickness;
        if (normalized >= 0.0 && normalized < options.wall_fraction)
            include(lower, profile[bin]);
        if (normalized <= 1.0 && normalized > 1.0 - options.wall_fraction)
            include(upper, profile[bin]);
        if (std::fabs(normalized - 0.5) <= 0.5 * options.core_fraction)
            include(core, profile[bin]);
    }
    RegionProfile wall = lower;
    add_profile(wall.profile, upper.profile);
    wall.volume += upper.volume;
    wall.bins += upper.bins;

    const auto density = [](const RegionProfile &region, double value) {
        return region.volume > 0.0 ? value / region.volume :
            std::numeric_limits<double>::quiet_NaN();
    };
    const auto total_mass = [](const RegionProfile &region) {
        return std::accumulate(region.profile.mass.begin(),
                               region.profile.mass.end(), 0.0);
    };
    const auto total_sites = [](const RegionProfile &region) {
        return std::accumulate(region.profile.functional_total.begin(),
                               region.profile.functional_total.end(), 0LL);
    };
    const auto reacted_sites = [](const RegionProfile &region) {
        return std::accumulate(region.profile.functional_reacted.begin(),
                               region.profile.functional_reacted.end(), 0LL);
    };
    const auto defect_count = [](const RegionProfile &region) {
        const BinCounts &topology = region.profile.topology;
        return topology.dangling_ends + topology.dangling_loops +
            topology.self_loops + topology.isolated_parents;
    };

    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "metric\tlower_wall\tupper_wall\tcombined_wall\tcore"
        << "\twall_over_core\tupper_minus_lower_relative\tunits\n";
    out << std::setprecision(12);
    const auto line = [&](const std::string &name, const std::string &units,
                          const auto &value) {
        const double lo = value(lower), hi = value(upper);
        const double wall_value = value(wall), core_value = value(core);
        const double wall_mean = 0.5 * (lo + hi);
        out << name << '\t' << lo << '\t' << hi << '\t' << wall_value
            << '\t' << core_value << '\t'
            << ratio_or_nan(wall_value, core_value) << '\t'
            << ratio_or_nan(hi - lo, wall_mean) << '\t' << units << '\n';
    };
    line("mass_density", "g/cm^3", [&](const RegionProfile &region) {
        return density(region, total_mass(region)) / kAvogadroAngstrom;
    });
    line("crosslink_bond_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, region.profile.topology.crosslink_bonds);
    });
    line("junction_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, region.profile.topology.junctions);
    });
    line("active_strand_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, region.profile.topology.active_strands);
    });
    line("defect_density", "A^-3", [&](const RegionProfile &region) {
        return density(region, defect_count(region));
    });
    line("local_conversion", "fraction", [&](const RegionProfile &region) {
        return ratio_or_nan(reacted_sites(region), total_sites(region));
    });
    line("active_contour_length_density", "A^-2",
         [&](const RegionProfile &region) {
        return density(region, region.profile.active_contour_length);
    });
    line("defect_contour_length_density", "A^-2",
         [&](const RegionProfile &region) {
        return density(region, region.profile.defect_contour_length);
    });
    line("segment_P2_z", "dimensionless", [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.segment_p2z_sum,
                           region.profile.contour_segments);
    });
    line("active_segment_P2_z", "dimensionless",
         [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_segment_p2z_sum,
                           region.profile.active_contour_segments);
    });
    line("active_strand_P2_z", "dimensionless",
         [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_p2z_sum,
                           region.profile.active_orientation_count);
    });
    line("active_Ree_parallel2_mean", "A^2",
         [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_ree_x2_sum +
                           region.profile.active_ree_y2_sum,
                           region.profile.active_orientation_count);
    });
    line("active_Ree_z2_mean", "A^2", [&](const RegionProfile &region) {
        return mean_or_nan(region.profile.active_ree_z2_sum,
                           region.profile.active_orientation_count);
    });
    if (z1_available) {
        line("z1_kink_density", "A^-3", [&](const RegionProfile &region) {
            return density(region, region.profile.z1_kinks);
        });
        line("z1_primitive_length_density", "A^-2",
             [&](const RegionProfile &region) {
            return density(region, region.profile.z1_primitive_length);
        });
        line("z1_primitive_P2_z", "dimensionless",
             [&](const RegionProfile &region) {
            return mean_or_nan(region.profile.z1_segment_p2z_sum,
                               region.profile.z1_segments);
        });
    }
}

struct DebyeWallerFrame {
    long long frame = 0;
    long long timestep = 0;
    double time_ps = 0.0;
    Vec3 drift;
    Vec3 global_msd;
    std::vector<Vec3> layer_msd;
};

double inverse_or_nan(double value) {
    return value > 0.0 ? 1.0 / value :
        std::numeric_limits<double>::quiet_NaN();
}

void write_debye_waller(
    const std::filesystem::path &global_path,
    const std::filesystem::path &layer_path,
    const Options &options, const ModelInfo &info, const DataFile &data) {
    DumpReader reader(options.debye_waller_trajectory_file);
    DumpFrame origin;
    if (!reader.next(origin, data.declared_atoms))
        throw std::runtime_error("Debye-Waller trajectory contains no frames");
    const int bins = std::max(1, static_cast<int>(
        std::ceil(origin.box.lz() / options.bin_width)));
    const double width = origin.box.lz() / bins;
    std::vector<int> origin_bin(origin.unwrapped.size(), -1);
    std::vector<long long> bin_counts(static_cast<std::size_t>(bins), 0);
    long long selected_beads = 0;
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        if (component_for_molecule(atom.molecule, info) != kStrand) continue;
        const int bin = bin_index(origin.unwrapped[static_cast<std::size_t>(id)].z,
                                  origin.box, bins, info.periodic_z());
        origin_bin[static_cast<std::size_t>(id)] = bin;
        ++bin_counts[static_cast<std::size_t>(bin)];
        ++selected_beads;
    }
    if (selected_beads == 0)
        throw std::runtime_error("Debye-Waller trajectory has no component-1 beads");

    std::vector<DebyeWallerFrame> frames;
    long long frame_index = 0;
    DumpFrame current = origin;
    while (true) {
        if (frame_index % options.frame_stride == 0) {
            DebyeWallerFrame sampled;
            sampled.frame = frame_index;
            sampled.timestep = current.timestep;
            sampled.time_ps =
                (current.timestep - origin.timestep) * info.timestep_fs * 1.0e-3;
            for (long long id = 1; id <= data.declared_atoms; ++id)
                sampled.drift += current.unwrapped[static_cast<std::size_t>(id)] -
                                 origin.unwrapped[static_cast<std::size_t>(id)];
            sampled.drift = (1.0 / data.declared_atoms) * sampled.drift;
            std::vector<Vec3> sums(static_cast<std::size_t>(bins));
            Vec3 global_sum;
            for (long long id = 1; id <= data.declared_atoms; ++id) {
                const int bin = origin_bin[static_cast<std::size_t>(id)];
                if (bin < 0) continue;
                const Vec3 displacement =
                    current.unwrapped[static_cast<std::size_t>(id)] -
                    origin.unwrapped[static_cast<std::size_t>(id)] - sampled.drift;
                Vec3 squared{displacement.x * displacement.x,
                             displacement.y * displacement.y,
                             displacement.z * displacement.z};
                sums[static_cast<std::size_t>(bin)] += squared;
                global_sum += squared;
            }
            sampled.global_msd = (1.0 / selected_beads) * global_sum;
            sampled.layer_msd.resize(static_cast<std::size_t>(bins));
            for (int bin = 0; bin < bins; ++bin) {
                const long long count = bin_counts[static_cast<std::size_t>(bin)];
                sampled.layer_msd[static_cast<std::size_t>(bin)] = count > 0
                    ? (1.0 / count) * sums[static_cast<std::size_t>(bin)]
                    : Vec3{std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()};
            }
            frames.push_back(std::move(sampled));
        }
        ++frame_index;
        if (!reader.next(current, data.declared_atoms)) break;
    }
    if (frames.size() < 3)
        throw std::runtime_error(
            "Debye-Waller trajectory needs at least three sampled frames");

    std::vector<double> logarithmic_slope(
        frames.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t index = 1; index + 1 < frames.size(); ++index) {
        const double first_time = frames[index - 1].time_ps;
        const double last_time = frames[index + 1].time_ps;
        const Vec3 first_msd = frames[index - 1].global_msd;
        const Vec3 last_msd = frames[index + 1].global_msd;
        const double first_total = first_msd.x + first_msd.y + first_msd.z;
        const double last_total = last_msd.x + last_msd.y + last_msd.z;
        if (first_time > 0.0 && last_time > first_time &&
            first_total > 0.0 && last_total > 0.0)
            logarithmic_slope[index] =
                std::log(last_total / first_total) /
                std::log(last_time / first_time);
    }

    std::size_t selected = 1;
    std::string selection_method;
    if (std::isfinite(options.debye_waller_time_ps)) {
        double nearest = std::numeric_limits<double>::infinity();
        for (std::size_t index = 1; index < frames.size(); ++index) {
            const double distance =
                std::fabs(frames[index].time_ps - options.debye_waller_time_ps);
            if (distance < nearest) {
                nearest = distance;
                selected = index;
            }
        }
        selection_method = "explicit_nearest_frame";
    } else {
        double minimum_slope = std::numeric_limits<double>::infinity();
        bool found = false;
        for (std::size_t index = 1; index + 1 < frames.size(); ++index) {
            if (frames[index].time_ps < options.debye_waller_search_start_ps ||
                frames[index].time_ps > options.debye_waller_search_end_ps ||
                !std::isfinite(logarithmic_slope[index])) continue;
            if (logarithmic_slope[index] < minimum_slope) {
                minimum_slope = logarithmic_slope[index];
                selected = index;
                found = true;
            }
        }
        if (found) {
            selection_method = "minimum_logarithmic_slope";
        } else {
            constexpr double kFallbackTimePs = 4.0;
            double nearest = std::numeric_limits<double>::infinity();
            for (std::size_t index = 1; index < frames.size(); ++index) {
                const double distance =
                    std::fabs(frames[index].time_ps - kFallbackTimePs);
                if (distance < nearest) {
                    nearest = distance;
                    selected = index;
                }
            }
            selection_method = "fallback_nearest_4ps";
        }
    }

    std::ofstream global(global_path);
    if (!global) throw std::runtime_error("cannot write " + global_path.string());
    global << "frame\ttimestep\ttime_ps\tstrand_beads"
        << "\tmsd_x_A2\tmsd_y_A2\tmsd_z_A2\tu2_xy_A2\tu2_3D_A2"
        << "\tlog_slope_u2_3D\tdrift_x_A\tdrift_y_A\tdrift_z_A\tselected\n";
    global << std::setprecision(12);
    for (std::size_t index = 0; index < frames.size(); ++index) {
        const DebyeWallerFrame &frame = frames[index];
        const double xy = frame.global_msd.x + frame.global_msd.y;
        global << frame.frame << '\t' << frame.timestep << '\t' << frame.time_ps
            << '\t' << selected_beads << '\t' << frame.global_msd.x << '\t'
            << frame.global_msd.y << '\t' << frame.global_msd.z << '\t'
            << xy << '\t' << xy + frame.global_msd.z << '\t'
            << logarithmic_slope[index] << '\t' << frame.drift.x << '\t'
            << frame.drift.y << '\t' << frame.drift.z << '\t'
            << (index == selected ? 1 : 0) << '\n';
    }

    const DebyeWallerFrame &dw = frames[selected];
    const double global_xy = dw.global_msd.x + dw.global_msd.y;
    const double global_total = global_xy + dw.global_msd.z;
    std::ofstream layer(layer_path);
    if (!layer) throw std::runtime_error("cannot write " + layer_path.string());
    layer << "bin\tzlo_origin_A\tzhi_origin_A\tstrand_beads"
        << "\tselected_time_ps\tselection_method\tlog_slope_u2_3D"
        << "\tu2_x_A2\tu2_y_A2\tu2_z_A2\tu2_xy_A2\tu2_3D_A2"
        << "\tu2_xy_over_global\tu2_3D_over_global"
        << "\tlocal_stiffness_xy_A-2\tlocal_stiffness_3D_A-2"
        << "\tstiffness_xy_over_global\tstiffness_3D_over_global\n";
    layer << std::setprecision(12);
    for (int bin = 0; bin < bins; ++bin) {
        const Vec3 value = dw.layer_msd[static_cast<std::size_t>(bin)];
        const double xy = value.x + value.y;
        const double total = xy + value.z;
        const double zlo = origin.box.zlo + bin * width;
        layer << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
            << bin_counts[static_cast<std::size_t>(bin)] << '\t' << dw.time_ps
            << '\t' << selection_method << '\t' << logarithmic_slope[selected]
            << '\t' << value.x << '\t' << value.y << '\t' << value.z << '\t'
            << xy << '\t' << total << '\t' << ratio_or_nan(xy, global_xy)
            << '\t' << ratio_or_nan(total, global_total) << '\t'
            << inverse_or_nan(xy) << '\t' << inverse_or_nan(total) << '\t'
            << ratio_or_nan(global_xy, xy) << '\t'
            << ratio_or_nan(global_total, total) << '\n';
    }
}

void write_layer_dynamics(
    const std::filesystem::path &msd_path,
    const std::filesystem::path &diffusion_path,
    const Options &options, const ModelInfo &info, const DataFile &data) {
    DumpReader reader(options.trajectory_file);
    DumpFrame origin;
    if (!reader.next(origin, data.declared_atoms))
        throw std::runtime_error("trajectory contains no frames");
    const int bins = std::max(1, static_cast<int>(
        std::ceil(origin.box.lz() / options.bin_width)));
    std::vector<int> origin_bin(origin.unwrapped.size(), -1);
    std::vector<long long> bin_counts(static_cast<std::size_t>(bins), 0);
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        if (component_for_molecule(atom.molecule, info) != kStrand) continue;
        const int bin = bin_index(origin.unwrapped[static_cast<std::size_t>(id)].z,
                                  origin.box, bins, info.periodic_z());
        origin_bin[static_cast<std::size_t>(id)] = bin;
        ++bin_counts[static_cast<std::size_t>(bin)];
    }
    std::ofstream out(msd_path);
    if (!out) throw std::runtime_error("cannot write " + msd_path.string());
    out << "frame\ttimestep\ttime_ns\tbin\tzlo_origin_A\tzhi_origin_A"
        << "\tstrand_beads\tmsd_x_A2\tmsd_y_A2\tmsd_z_A2"
        << "\tmsd_parallel_A2\tmsd_total_A2\tdrift_x_A\tdrift_y_A\tdrift_z_A\n";
    const double width = origin.box.lz() / bins;
    std::vector<double> sampled_times;
    std::vector<std::vector<double>> parallel_msd(static_cast<std::size_t>(bins));
    std::vector<std::vector<double>> total_msd(static_cast<std::size_t>(bins));
    long long frame_index = 0;
    DumpFrame current = origin;
    while (true) {
        if (frame_index % options.frame_stride == 0) {
            Vec3 drift;
            for (long long id = 1; id <= data.declared_atoms; ++id)
                drift += current.unwrapped[static_cast<std::size_t>(id)] -
                         origin.unwrapped[static_cast<std::size_t>(id)];
            drift = (1.0 / data.declared_atoms) * drift;
            std::vector<Vec3> sum(static_cast<std::size_t>(bins));
            for (long long id = 1; id <= data.declared_atoms; ++id) {
                const int bin = origin_bin[static_cast<std::size_t>(id)];
                if (bin < 0) continue;
                const Vec3 displacement =
                    current.unwrapped[static_cast<std::size_t>(id)] -
                    origin.unwrapped[static_cast<std::size_t>(id)] - drift;
                sum[static_cast<std::size_t>(bin)].x += displacement.x * displacement.x;
                sum[static_cast<std::size_t>(bin)].y += displacement.y * displacement.y;
                sum[static_cast<std::size_t>(bin)].z += displacement.z * displacement.z;
            }
            const double time_ns =
                (current.timestep - origin.timestep) * info.timestep_fs * 1.0e-6;
            sampled_times.push_back(time_ns);
            for (int bin = 0; bin < bins; ++bin) {
                const long long count = bin_counts[static_cast<std::size_t>(bin)];
                const Vec3 msd = count > 0
                    ? (1.0 / count) * sum[static_cast<std::size_t>(bin)] :
                      Vec3{std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::quiet_NaN()};
                const double zlo = origin.box.zlo + bin * width;
                const double parallel = msd.x + msd.y;
                const double total = parallel + msd.z;
                parallel_msd[static_cast<std::size_t>(bin)].push_back(parallel);
                total_msd[static_cast<std::size_t>(bin)].push_back(total);
                out << frame_index << '\t' << current.timestep << '\t'
                    << std::setprecision(12) << time_ns << '\t' << bin + 1 << '\t'
                    << zlo << '\t' << zlo + width << '\t' << count << '\t'
                    << msd.x << '\t' << msd.y << '\t' << msd.z << '\t'
                    << parallel << '\t' << total << '\t'
                    << drift.x << '\t' << drift.y << '\t' << drift.z << '\n';
            }
        }
        ++frame_index;
        if (!reader.next(current, data.declared_atoms)) break;
    }

    std::ofstream diffusion(diffusion_path);
    if (!diffusion)
        throw std::runtime_error("cannot write " + diffusion_path.string());
    diffusion
        << "bin\tzlo_origin_A\tzhi_origin_A\tstrand_beads"
        << "\tfit_start_ns\tfit_end_ns"
        << "\tfit_points_xy\tslope_xy_A2_per_ns\tintercept_xy_A2\tR2_xy"
        << "\tD_xy_A2_per_ns\tD_xy_cm2_per_s"
        << "\tfit_points_3D\tslope_3D_A2_per_ns\tintercept_3D_A2\tR2_3D"
        << "\tD_3D_A2_per_ns\tD_3D_cm2_per_s"
        << "\trecommended_dimension\trecommended_D_A2_per_ns"
        << "\trecommended_D_cm2_per_s\trecommended_R2\n";
    const double end_time = sampled_times.empty() ? 0.0 : sampled_times.back();
    const double start_time = options.diffusion_fit_start_fraction * end_time;
    constexpr double kAngstrom2PerNsToCm2PerS = 1.0e-7;
    for (int bin = 0; bin < bins; ++bin) {
        const LinearFit xy = linear_fit(
            sampled_times, parallel_msd[static_cast<std::size_t>(bin)], start_time);
        const LinearFit three_d = linear_fit(
            sampled_times, total_msd[static_cast<std::size_t>(bin)], start_time);
        const double d_xy = xy.slope / 4.0;
        const double d_3d = three_d.slope / 6.0;
        const bool recommend_xy = info.geometry == "film";
        const double recommended_d = recommend_xy ? d_xy : d_3d;
        const double recommended_r2 = recommend_xy ? xy.r_squared : three_d.r_squared;
        const double zlo = origin.box.zlo + bin * width;
        diffusion << bin + 1 << '\t' << std::setprecision(12)
            << zlo << '\t' << zlo + width << '\t'
            << bin_counts[static_cast<std::size_t>(bin)] << '\t'
            << start_time << '\t' << end_time << '\t'
            << xy.points << '\t' << xy.slope << '\t' << xy.intercept << '\t'
            << xy.r_squared << '\t' << d_xy << '\t'
            << d_xy * kAngstrom2PerNsToCm2PerS << '\t'
            << three_d.points << '\t' << three_d.slope << '\t'
            << three_d.intercept << '\t' << three_d.r_squared << '\t'
            << d_3d << '\t' << d_3d * kAngstrom2PerNsToCm2PerS << '\t'
            << (recommend_xy ? "2D_xy" : "3D") << '\t'
            << recommended_d << '\t'
            << recommended_d * kAngstrom2PerNsToCm2PerS << '\t'
            << recommended_r2 << '\n';
    }
}

struct TimeOriginFrame {
    long long timestep = 0;
    Vec3 system_center;
    std::vector<Vec3> strand_positions;
    std::vector<int> strand_bins;
    std::vector<long long> bin_counts;
};

struct TimeAveragedLag {
    explicit TimeAveragedLag(int bins = 0)
        : squared_sum(static_cast<std::size_t>(bins)),
          observations(static_cast<std::size_t>(bins), 0) {}

    long long time_origins = 0;
    Vec3 drift_sum;
    std::vector<Vec3> squared_sum;
    std::vector<long long> observations;
};

Vec3 system_center(const DumpFrame &frame, long long atom_count) {
    Vec3 center;
    for (long long id = 1; id <= atom_count; ++id)
        center += frame.unwrapped[static_cast<std::size_t>(id)];
    return (1.0 / atom_count) * center;
}

void write_time_averaged_debye_waller(
    const std::filesystem::path &global_path,
    const std::filesystem::path &layer_path,
    const Options &options, const ModelInfo &info, const DataFile &data) {
    DumpReader reader(options.debye_waller_trajectory_file);
    DumpFrame reference;
    if (!reader.next(reference, data.declared_atoms))
        throw std::runtime_error("Debye-Waller trajectory contains no frames");
    const int bins = std::max(1, static_cast<int>(
        std::ceil(reference.box.lz() / options.bin_width)));
    const double width = reference.box.lz() / bins;

    std::vector<long long> strand_ids;
    strand_ids.reserve(static_cast<std::size_t>(data.declared_atoms));
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        if (component_for_molecule(atom.molecule, info) == kStrand)
            strand_ids.push_back(id);
    }
    if (strand_ids.empty())
        throw std::runtime_error(
            "Debye-Waller trajectory has no component-1 beads");

    std::vector<TimeOriginFrame> origins;
    origins.reserve(static_cast<std::size_t>(options.dw_time_origin_count));
    std::map<long long, TimeAveragedLag> lag_accumulators;
    long long frame_index = 0;
    long long sampled_frame_index = 0;
    DumpFrame current = reference;
    while (true) {
        if (frame_index % options.frame_stride == 0) {
            const Vec3 current_center =
                system_center(current, data.declared_atoms);
            if (sampled_frame_index % options.dw_time_origin_stride == 0 &&
                static_cast<long long>(origins.size()) <
                    options.dw_time_origin_count) {
                TimeOriginFrame origin;
                origin.timestep = current.timestep;
                origin.system_center = current_center;
                origin.strand_positions.reserve(strand_ids.size());
                origin.strand_bins.reserve(strand_ids.size());
                origin.bin_counts.assign(static_cast<std::size_t>(bins), 0);
                for (const long long id : strand_ids) {
                    const Vec3 position =
                        current.unwrapped[static_cast<std::size_t>(id)];
                    const int bin = bin_index(
                        position.z, reference.box, bins, info.periodic_z());
                    origin.strand_positions.push_back(position);
                    origin.strand_bins.push_back(bin);
                    ++origin.bin_counts[static_cast<std::size_t>(bin)];
                }
                origins.push_back(std::move(origin));
            }

            for (const TimeOriginFrame &origin : origins) {
                const long long lag_steps = current.timestep - origin.timestep;
                if (lag_steps < 0)
                    throw std::runtime_error(
                        "Debye-Waller timesteps are not monotonically increasing");
                auto inserted = lag_accumulators.try_emplace(lag_steps, bins);
                TimeAveragedLag &lag = inserted.first->second;
                ++lag.time_origins;
                const Vec3 drift = current_center - origin.system_center;
                lag.drift_sum += drift;
                for (int bin = 0; bin < bins; ++bin)
                    lag.observations[static_cast<std::size_t>(bin)] +=
                        origin.bin_counts[static_cast<std::size_t>(bin)];
                for (std::size_t index = 0; index < strand_ids.size(); ++index) {
                    const long long id = strand_ids[index];
                    const int bin = origin.strand_bins[index];
                    const Vec3 displacement =
                        current.unwrapped[static_cast<std::size_t>(id)] -
                        origin.strand_positions[index] - drift;
                    Vec3 &sum = lag.squared_sum[static_cast<std::size_t>(bin)];
                    sum.x += displacement.x * displacement.x;
                    sum.y += displacement.y * displacement.y;
                    sum.z += displacement.z * displacement.z;
                }
            }
            ++sampled_frame_index;
        }
        ++frame_index;
        if (!reader.next(current, data.declared_atoms)) break;
    }
    if (origins.size() < 2)
        throw std::runtime_error(
            "time-averaged Debye-Waller analysis needs at least two origins; "
            "reduce --dw-time-origin-stride or increase "
            "--dw-time-origin-count");
    if (lag_accumulators.size() < 3)
        throw std::runtime_error(
            "time-averaged Debye-Waller analysis needs at least three lags");

    std::vector<long long> lag_steps;
    std::vector<const TimeAveragedLag *> lags;
    std::vector<double> time_ps;
    std::vector<Vec3> global_msd;
    std::vector<long long> global_observations;
    for (const auto &entry : lag_accumulators) {
        lag_steps.push_back(entry.first);
        lags.push_back(&entry.second);
        time_ps.push_back(entry.first * info.timestep_fs * 1.0e-3);
        Vec3 sum;
        long long observations = 0;
        for (int bin = 0; bin < bins; ++bin) {
            sum += entry.second.squared_sum[static_cast<std::size_t>(bin)];
            observations +=
                entry.second.observations[static_cast<std::size_t>(bin)];
        }
        global_observations.push_back(observations);
        global_msd.push_back(observations > 0
            ? (1.0 / observations) * sum
            : Vec3{std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::quiet_NaN()});
    }

    std::vector<double> logarithmic_slope(
        lags.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t index = 1; index + 1 < lags.size(); ++index) {
        const double first_total = global_msd[index - 1].x +
            global_msd[index - 1].y + global_msd[index - 1].z;
        const double last_total = global_msd[index + 1].x +
            global_msd[index + 1].y + global_msd[index + 1].z;
        if (time_ps[index - 1] > 0.0 &&
            time_ps[index + 1] > time_ps[index - 1] &&
            first_total > 0.0 && last_total > 0.0)
            logarithmic_slope[index] =
                std::log(last_total / first_total) /
                std::log(time_ps[index + 1] / time_ps[index - 1]);
    }

    std::size_t selected = 1;
    std::string selection_method;
    if (std::isfinite(options.debye_waller_time_ps)) {
        double nearest = std::numeric_limits<double>::infinity();
        for (std::size_t index = 1; index < lags.size(); ++index) {
            const double distance =
                std::fabs(time_ps[index] - options.debye_waller_time_ps);
            if (distance < nearest) {
                nearest = distance;
                selected = index;
            }
        }
        selection_method = "explicit_nearest_lag";
    } else {
        double minimum_slope = std::numeric_limits<double>::infinity();
        bool found = false;
        for (std::size_t index = 1; index + 1 < lags.size(); ++index) {
            if (time_ps[index] < options.debye_waller_search_start_ps ||
                time_ps[index] > options.debye_waller_search_end_ps ||
                !std::isfinite(logarithmic_slope[index])) continue;
            if (logarithmic_slope[index] < minimum_slope) {
                minimum_slope = logarithmic_slope[index];
                selected = index;
                found = true;
            }
        }
        if (found) {
            selection_method = "minimum_logarithmic_slope";
        } else {
            constexpr double kFallbackTimePs = 4.0;
            double nearest = std::numeric_limits<double>::infinity();
            for (std::size_t index = 1; index < lags.size(); ++index) {
                const double distance =
                    std::fabs(time_ps[index] - kFallbackTimePs);
                if (distance < nearest) {
                    nearest = distance;
                    selected = index;
                }
            }
            selection_method = "fallback_nearest_4ps";
        }
    }

    std::ofstream global(global_path);
    if (!global) throw std::runtime_error("cannot write " + global_path.string());
    global << "lag_index\tlag_steps\ttime_ps\ttime_origins"
        << "\tstrand_bead_observations"
        << "\tmsd_x_A2\tmsd_y_A2\tmsd_z_A2\tu2_xy_A2\tu2_3D_A2"
        << "\tlog_slope_u2_3D\tmean_drift_x_A\tmean_drift_y_A"
        << "\tmean_drift_z_A\tselected\n";
    global << std::setprecision(12);
    for (std::size_t index = 0; index < lags.size(); ++index) {
        const Vec3 value = global_msd[index];
        const double xy = value.x + value.y;
        const Vec3 mean_drift =
            (1.0 / lags[index]->time_origins) * lags[index]->drift_sum;
        global << index << '\t' << lag_steps[index] << '\t' << time_ps[index]
            << '\t' << lags[index]->time_origins << '\t'
            << global_observations[index] << '\t' << value.x << '\t'
            << value.y << '\t' << value.z << '\t' << xy << '\t'
            << xy + value.z << '\t' << logarithmic_slope[index] << '\t'
            << mean_drift.x << '\t' << mean_drift.y << '\t'
            << mean_drift.z << '\t' << (index == selected ? 1 : 0) << '\n';
    }

    const TimeAveragedLag &dw = *lags[selected];
    const double global_xy =
        global_msd[selected].x + global_msd[selected].y;
    const double global_total = global_xy + global_msd[selected].z;
    std::ofstream layer(layer_path);
    if (!layer) throw std::runtime_error("cannot write " + layer_path.string());
    layer << "bin\tzlo_origin_A\tzhi_origin_A\ttime_origins"
        << "\tstrand_bead_observations\tmean_strand_beads_per_origin"
        << "\tselected_time_ps\tselection_method\tlog_slope_u2_3D"
        << "\tu2_x_A2\tu2_y_A2\tu2_z_A2\tu2_xy_A2\tu2_3D_A2"
        << "\tu2_xy_over_global\tu2_3D_over_global"
        << "\tlocal_stiffness_xy_A-2\tlocal_stiffness_3D_A-2"
        << "\tstiffness_xy_over_global\tstiffness_3D_over_global\n";
    layer << std::setprecision(12);
    for (int bin = 0; bin < bins; ++bin) {
        const long long observations =
            dw.observations[static_cast<std::size_t>(bin)];
        const Vec3 value = observations > 0
            ? (1.0 / observations) *
                dw.squared_sum[static_cast<std::size_t>(bin)]
            : Vec3{std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::quiet_NaN(),
                   std::numeric_limits<double>::quiet_NaN()};
        const double xy = value.x + value.y;
        const double total = xy + value.z;
        const double zlo = reference.box.zlo + bin * width;
        layer << bin + 1 << '\t' << zlo << '\t' << zlo + width << '\t'
            << dw.time_origins << '\t' << observations << '\t'
            << static_cast<double>(observations) / dw.time_origins << '\t'
            << time_ps[selected] << '\t' << selection_method << '\t'
            << logarithmic_slope[selected] << '\t' << value.x << '\t'
            << value.y << '\t' << value.z << '\t' << xy << '\t' << total
            << '\t' << ratio_or_nan(xy, global_xy) << '\t'
            << ratio_or_nan(total, global_total) << '\t'
            << inverse_or_nan(xy) << '\t' << inverse_or_nan(total) << '\t'
            << ratio_or_nan(global_xy, xy) << '\t'
            << ratio_or_nan(global_total, total) << '\n';
    }
}

void write_time_averaged_layer_dynamics(
    const std::filesystem::path &msd_path,
    const std::filesystem::path &diffusion_path,
    const Options &options, const ModelInfo &info, const DataFile &data) {
    DumpReader reader(options.trajectory_file);
    DumpFrame reference;
    if (!reader.next(reference, data.declared_atoms))
        throw std::runtime_error("trajectory contains no frames");
    const int bins = std::max(1, static_cast<int>(
        std::ceil(reference.box.lz() / options.bin_width)));
    const double width = reference.box.lz() / bins;

    std::vector<long long> strand_ids;
    strand_ids.reserve(static_cast<std::size_t>(data.declared_atoms));
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        if (component_for_molecule(atom.molecule, info) == kStrand)
            strand_ids.push_back(id);
    }
    if (strand_ids.empty())
        throw std::runtime_error("trajectory has no component-1 beads");

    std::vector<TimeOriginFrame> origins;
    origins.reserve(static_cast<std::size_t>(options.time_origin_count));
    std::map<long long, TimeAveragedLag> lag_accumulators;
    long long frame_index = 0;
    long long sampled_frame_index = 0;
    DumpFrame current = reference;
    while (true) {
        if (frame_index % options.frame_stride == 0) {
            const Vec3 current_center =
                system_center(current, data.declared_atoms);
            if (sampled_frame_index % options.time_origin_stride == 0 &&
                static_cast<long long>(origins.size()) <
                    options.time_origin_count) {
                TimeOriginFrame origin;
                origin.timestep = current.timestep;
                origin.system_center = current_center;
                origin.strand_positions.reserve(strand_ids.size());
                origin.strand_bins.reserve(strand_ids.size());
                origin.bin_counts.assign(static_cast<std::size_t>(bins), 0);
                for (const long long id : strand_ids) {
                    const Vec3 position =
                        current.unwrapped[static_cast<std::size_t>(id)];
                    const int bin = bin_index(
                        position.z, reference.box, bins, info.periodic_z());
                    origin.strand_positions.push_back(position);
                    origin.strand_bins.push_back(bin);
                    ++origin.bin_counts[static_cast<std::size_t>(bin)];
                }
                origins.push_back(std::move(origin));
            }

            for (const TimeOriginFrame &origin : origins) {
                const long long lag_steps = current.timestep - origin.timestep;
                if (lag_steps < 0)
                    throw std::runtime_error(
                        "trajectory timesteps are not monotonically increasing");
                auto inserted = lag_accumulators.try_emplace(lag_steps, bins);
                TimeAveragedLag &lag = inserted.first->second;
                ++lag.time_origins;
                const Vec3 drift = current_center - origin.system_center;
                lag.drift_sum += drift;
                for (int bin = 0; bin < bins; ++bin)
                    lag.observations[static_cast<std::size_t>(bin)] +=
                        origin.bin_counts[static_cast<std::size_t>(bin)];
                for (std::size_t index = 0; index < strand_ids.size(); ++index) {
                    const long long id = strand_ids[index];
                    const int bin = origin.strand_bins[index];
                    const Vec3 displacement =
                        current.unwrapped[static_cast<std::size_t>(id)] -
                        origin.strand_positions[index] - drift;
                    Vec3 &sum = lag.squared_sum[static_cast<std::size_t>(bin)];
                    sum.x += displacement.x * displacement.x;
                    sum.y += displacement.y * displacement.y;
                    sum.z += displacement.z * displacement.z;
                }
            }
            ++sampled_frame_index;
        }
        ++frame_index;
        if (!reader.next(current, data.declared_atoms)) break;
    }
    if (origins.size() < 2)
        throw std::runtime_error(
            "time-averaged layer MSD needs at least two selected time origins; "
            "reduce --time-origin-stride or increase --time-origin-count");

    std::ofstream out(msd_path);
    if (!out) throw std::runtime_error("cannot write " + msd_path.string());
    out << "lag_index\tlag_steps\ttime_ns\tbin\tzlo_origin_A\tzhi_origin_A"
        << "\ttime_origins\tstrand_bead_observations"
        << "\tmean_strand_beads_per_origin"
        << "\tmsd_x_A2\tmsd_y_A2\tmsd_z_A2"
        << "\tmsd_parallel_A2\tmsd_total_A2"
        << "\tmean_drift_x_A\tmean_drift_y_A\tmean_drift_z_A\n";

    std::vector<double> sampled_times;
    std::vector<long long> time_origins_by_lag;
    std::vector<std::vector<double>> parallel_msd(
        static_cast<std::size_t>(bins));
    std::vector<std::vector<double>> total_msd(
        static_cast<std::size_t>(bins));
    long long lag_index = 0;
    for (const auto &entry : lag_accumulators) {
        const long long lag_steps = entry.first;
        const TimeAveragedLag &lag = entry.second;
        const double time_ns =
            lag_steps * info.timestep_fs * 1.0e-6;
        sampled_times.push_back(time_ns);
        time_origins_by_lag.push_back(lag.time_origins);
        const Vec3 mean_drift =
            (1.0 / lag.time_origins) * lag.drift_sum;
        for (int bin = 0; bin < bins; ++bin) {
            const long long observations =
                lag.observations[static_cast<std::size_t>(bin)];
            const Vec3 msd = observations > 0
                ? (1.0 / observations) *
                    lag.squared_sum[static_cast<std::size_t>(bin)]
                : Vec3{std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN()};
            const double parallel = msd.x + msd.y;
            const double total = parallel + msd.z;
            parallel_msd[static_cast<std::size_t>(bin)].push_back(parallel);
            total_msd[static_cast<std::size_t>(bin)].push_back(total);
            const double zlo = reference.box.zlo + bin * width;
            out << lag_index << '\t' << lag_steps << '\t'
                << std::setprecision(12) << time_ns << '\t' << bin + 1 << '\t'
                << zlo << '\t' << zlo + width << '\t' << lag.time_origins
                << '\t' << observations << '\t'
                << static_cast<double>(observations) / lag.time_origins << '\t'
                << msd.x << '\t' << msd.y << '\t' << msd.z << '\t'
                << parallel << '\t' << total << '\t'
                << mean_drift.x << '\t' << mean_drift.y << '\t'
                << mean_drift.z << '\n';
        }
        ++lag_index;
    }

    std::ofstream diffusion(diffusion_path);
    if (!diffusion)
        throw std::runtime_error("cannot write " + diffusion_path.string());
    diffusion
        << "bin\tzlo_origin_A\tzhi_origin_A\tmean_strand_beads_per_origin"
        << "\tselected_time_origins\tminimum_time_origins_in_fit"
        << "\tfit_start_ns\tfit_end_ns"
        << "\tfit_points_xy\tslope_xy_A2_per_ns\tintercept_xy_A2\tR2_xy"
        << "\tD_xy_A2_per_ns\tD_xy_cm2_per_s"
        << "\tfit_points_3D\tslope_3D_A2_per_ns\tintercept_3D_A2\tR2_3D"
        << "\tD_3D_A2_per_ns\tD_3D_cm2_per_s"
        << "\trecommended_dimension\trecommended_D_A2_per_ns"
        << "\trecommended_D_cm2_per_s\trecommended_R2\n";
    const long long minimum_fit_origins =
        (static_cast<long long>(origins.size()) + 1) / 2;
    double end_time = 0.0;
    for (std::size_t index = 0; index < sampled_times.size(); ++index)
        if (time_origins_by_lag[index] >= minimum_fit_origins)
            end_time = sampled_times[index];
    const double start_time = options.diffusion_fit_start_fraction * end_time;
    constexpr double kAngstrom2PerNsToCm2PerS = 1.0e-7;
    const TimeAveragedLag &zero_lag = lag_accumulators.at(0);
    for (int bin = 0; bin < bins; ++bin) {
        const LinearFit xy = linear_fit(
            sampled_times, parallel_msd[static_cast<std::size_t>(bin)],
            start_time, end_time);
        const LinearFit three_d = linear_fit(
            sampled_times, total_msd[static_cast<std::size_t>(bin)],
            start_time, end_time);
        const double d_xy = xy.slope / 4.0;
        const double d_3d = three_d.slope / 6.0;
        const bool recommend_xy = info.geometry == "film";
        const double recommended_d = recommend_xy ? d_xy : d_3d;
        const double recommended_r2 =
            recommend_xy ? xy.r_squared : three_d.r_squared;
        const double zlo = reference.box.zlo + bin * width;
        const double mean_beads =
            static_cast<double>(
                zero_lag.observations[static_cast<std::size_t>(bin)]) /
            zero_lag.time_origins;
        diffusion << bin + 1 << '\t' << std::setprecision(12)
            << zlo << '\t' << zlo + width << '\t' << mean_beads << '\t'
            << origins.size() << '\t' << minimum_fit_origins << '\t'
            << start_time << '\t' << end_time << '\t'
            << xy.points << '\t' << xy.slope << '\t' << xy.intercept << '\t'
            << xy.r_squared << '\t' << d_xy << '\t'
            << d_xy * kAngstrom2PerNsToCm2PerS << '\t'
            << three_d.points << '\t' << three_d.slope << '\t'
            << three_d.intercept << '\t' << three_d.r_squared << '\t'
            << d_3d << '\t' << d_3d * kAngstrom2PerNsToCm2PerS << '\t'
            << (recommend_xy ? "2D_xy" : "3D") << '\t'
            << recommended_d << '\t'
            << recommended_d * kAngstrom2PerNsToCm2PerS << '\t'
            << recommended_r2 << '\n';
    }
}

void write_report(
    const std::filesystem::path &path, const Options &options,
    const ModelInfo &info, const DataFile &data, int bins,
    const std::string &z1_sp_file) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    const SlidingGrid sliding_grid = make_sliding_grid(data.box, options);
    out << "PDMS network distribution and dynamics report\n"
        << "case: " << info.case_name << "\n"
        << "geometry: " << info.geometry << "\n"
        << "z bins: " << bins << "\n"
        << "realized bin width: " << data.box.lz() / bins << " A\n"
        << "static sliding profile requested window / step: "
        << options.profile_window_width << " / " << options.profile_step
        << " A\n"
        << "static sliding profile realized window / step: "
        << sliding_grid.window_width_A << " / " << sliding_grid.step_A
        << " A\n"
        << "static sliding windows overlap and are not independent samples\n"
        << "static markers: reaction-bond midpoints, reduced-strand midpoints,"
        << " free dangling ends, loop markers, and parent centers\n"
        << "densities use the full lateral box area and physical bin volume\n"
        << "local conversion: reacted functional sites divided by declared"
        << " functional sites at their final z positions\n"
        << "crosslink bond position: periodicity-corrected final bond midpoint\n"
        << "contour profile: reduced-strand segment length assigned by segment midpoint\n"
        << "defect contour: every reduced strand not classified active\n"
        << "segment orientation: P2_alpha=(3*u_alpha^2-1)/2\n"
        << "strand conformation: assigned by contour midpoint\n";
    if (info.geometry == "film") {
        double thickness = info.film_thickness_angstrom;
        if (!(thickness > 0.0))
            thickness = data.box.lz() -
                2.0 * info.film_wall_cutoff_per_side_angstrom;
        out << "film box Lz: " << data.box.lz() << " A\n"
            << "film wall cutoff per side: "
            << info.film_wall_cutoff_per_side_angstrom << " A\n"
            << "nominal wall-free material thickness: " << thickness << " A\n"
            << "material coordinate: (z-zlo-wall_cutoff)/material_thickness\n"
            << "summary wall fraction per side: " << options.wall_fraction << "\n"
            << "summary centered core fraction: " << options.core_fraction << "\n";
    } else {
        out << "bulk profile guidance: z is periodic and has no physical wall;"
            << " folded and boundary/core outputs are uniformity controls only\n";
    }
    if (z1_sp_file.empty()) {
        out << "Z1+ profile: unavailable; z1 density/orientation columns are NaN\n";
    } else {
        out << "Z1+ primitive-path result: " << z1_sp_file << "\n"
            << "Z1+ kink definition: primitive-path points with nonzero kink flag\n"
            << "Z1+ scaling: none; result box lengths must match the snapshot\n";
    }
    if (options.debye_waller_trajectory_file.empty()) {
        out << "Debye-Waller analysis: not requested\n";
    } else {
        out << "Debye-Waller trajectory: "
            << options.debye_waller_trajectory_file << "\n"
            << "Debye-Waller selection: component-1 layer MSD with whole-system"
            << " COM drift removed\n";
        if (std::isfinite(options.debye_waller_time_ps))
            out << "Debye-Waller requested time: "
                << options.debye_waller_time_ps << " ps (nearest frame)\n";
        else
            out << "Debye-Waller time: minimum global 3D MSD logarithmic slope"
                << " from " << options.debye_waller_search_start_ps << " to "
                << options.debye_waller_search_end_ps << " ps\n";
        if (options.time_averaged_dw) {
            out << "time-averaged Debye-Waller analysis: enabled\n"
                << "Debye-Waller time-origin layer assignment: component-1 beads"
                << " are reassigned by z at each selected origin\n"
                << "Debye-Waller time-origin selection: at most "
                << options.dw_time_origin_count << " origins, every "
                << options.dw_time_origin_stride << " sampled frames\n";
        } else {
            out << "time-averaged Debye-Waller analysis: not requested\n";
        }
    }
    if (options.trajectory_file.empty()) {
        out << "layer dynamics: not requested\n";
    } else {
        out << "trajectory: " << options.trajectory_file << "\n"
            << "layer dynamics selection: all component-1 beads, grouped by first-frame z\n"
            << "MSD convention: displacement from first frame with whole-system COM drift removed\n"
            << "trajectory-layer bins: derived from first-frame box and requested bin width\n"
            << "frame stride: " << options.frame_stride << "\n"
            << "diffusion fit start fraction: "
            << options.diffusion_fit_start_fraction << " of recorded duration\n"
            << "D_xy definition: slope(MSD_x+MSD_y)/4\n"
            << "D_3D definition: slope(MSD_total)/6\n"
            << "diffusion slopes remain signed; negative values flag unresolved diffusion or noise\n";
        if (options.time_averaged_msd) {
            out << "time-averaged layer MSD: enabled\n"
                << "time-origin layer assignment: component-1 beads are reassigned by z at each selected origin\n"
                << "time-origin selection: at most "
                << options.time_origin_count << " origins, every "
                << options.time_origin_stride << " sampled frames\n"
                << "time-averaged lag grouping: exact timestep difference\n"
                << "time-averaged drift correction: whole-system COM displacement for each origin pair\n"
                << "time-averaged diffusion fit end: largest lag sampled by at least half of selected origins\n";
        } else {
            out << "time-averaged layer MSD: not requested\n";
        }
    }
    if (info.geometry == "film")
        out << "film guidance: compare x/y or parallel MSD; interpret z using the"
            << " trajectory boundary condition (wall-confined or free-surface)\n";
}

void print_help(const char *program) {
    std::cout
        << "Usage: " << program << " <case>.npt_eq <case>.info [options]\n\n"
        << "Options:\n"
        << "  --dw-trajectory FILE optional high-frequency Debye-Waller dump\n"
        << "  --trajectory FILE    optional dump.msd.lammpstrj for layer MSD\n"
        << "  --z1-sp FILE         override the default Z1+SP.dat path\n"
        << "  --no-z1              disable default Z1+SP.dat auto-detection\n"
        << "  --bin-width X        target z-bin width in A (default 5)\n"
        << "  --profile-window X   sliding static-profile window in A (default 5)\n"
        << "  --profile-step X     sliding static-profile center spacing in A (default 1)\n"
        << "  --wall-fraction X    material fraction at each wall for summary (default 0.20)\n"
        << "  --core-fraction X    centered material fraction for summary (default 0.20)\n"
        << "  --frame-stride N     analyze every Nth trajectory frame (default 1)\n"
        << "  --time-averaged-msd  also calculate layer MSD over sampled time origins\n"
        << "  --time-origin-stride N\n"
        << "                      spacing between time origins in sampled frames (default 100)\n"
        << "  --time-origin-count N\n"
        << "                      maximum selected time origins (default 10)\n"
        << "  --time-averaged-dw   also calculate Debye-Waller MSD over sampled origins\n"
        << "  --dw-time-origin-stride N\n"
        << "                      DW origin spacing in sampled frames (default 100)\n"
        << "  --dw-time-origin-count N\n"
        << "                      maximum selected DW origins (default 10)\n"
        << "  --diffusion-fit-start-fraction X\n"
        << "                      fit D over final 1-X fraction (default 0.50)\n"
        << "  --dw-time-ps X       use the nearest DW frame to X ps\n"
        << "  --dw-search-start-ps X\n"
        << "                      automatic DW search start (default 0.50 ps)\n"
        << "  --dw-search-end-ps X automatic DW search end (default 20 ps)\n"
        << "  --output-dir PATH    output directory (default analysis_<case>)\n"
        << "  --help               show this help\n";
}

Options parse_options(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        print_help(argv[0]); std::exit(0);
    }
    if (argc < 3) throw std::runtime_error("expected data and info files");
    Options options;
    options.data_file = argv[1];
    options.info_file = argv[2];
    for (int index = 3; index < argc; ++index) {
        const std::string option = argv[index];
        const auto value = [&]() {
            if (++index >= argc) throw std::runtime_error("missing value for " + option);
            return std::string(argv[index]);
        };
        if (option == "--dw-trajectory")
            options.debye_waller_trajectory_file = value();
        else if (option == "--trajectory") options.trajectory_file = value();
        else if (option == "--z1-sp" || option == "--z1-results")
            options.z1_sp_file = value();
        else if (option == "--no-z1") options.disable_z1 = true;
        else if (option == "--bin-width") options.bin_width = std::stod(value());
        else if (option == "--profile-window")
            options.profile_window_width = std::stod(value());
        else if (option == "--profile-step")
            options.profile_step = std::stod(value());
        else if (option == "--wall-fraction")
            options.wall_fraction = std::stod(value());
        else if (option == "--core-fraction")
            options.core_fraction = std::stod(value());
        else if (option == "--frame-stride") options.frame_stride = std::stoll(value());
        else if (option == "--time-averaged-msd")
            options.time_averaged_msd = true;
        else if (option == "--time-origin-stride")
            options.time_origin_stride = std::stoll(value());
        else if (option == "--time-origin-count")
            options.time_origin_count = std::stoll(value());
        else if (option == "--time-averaged-dw")
            options.time_averaged_dw = true;
        else if (option == "--dw-time-origin-stride")
            options.dw_time_origin_stride = std::stoll(value());
        else if (option == "--dw-time-origin-count")
            options.dw_time_origin_count = std::stoll(value());
        else if (option == "--diffusion-fit-start-fraction")
            options.diffusion_fit_start_fraction = std::stod(value());
        else if (option == "--dw-time-ps")
            options.debye_waller_time_ps = std::stod(value());
        else if (option == "--dw-search-start-ps")
            options.debye_waller_search_start_ps = std::stod(value());
        else if (option == "--dw-search-end-ps")
            options.debye_waller_search_end_ps = std::stod(value());
        else if (option == "--output-dir") options.output_directory = value();
        else if (option == "--help") { print_help(argv[0]); std::exit(0); }
        else throw std::runtime_error("unknown option: " + option);
    }
    if (options.bin_width <= 0.0) throw std::runtime_error("bin width must be positive");
    if (!(options.profile_window_width > 0.0))
        throw std::runtime_error("profile window must be positive");
    if (!(options.profile_step > 0.0))
        throw std::runtime_error("profile step must be positive");
    if (!(options.wall_fraction > 0.0 && options.wall_fraction < 0.5))
        throw std::runtime_error("wall fraction must be between 0 and 0.5");
    if (!(options.core_fraction > 0.0 && options.core_fraction <= 1.0))
        throw std::runtime_error("core fraction must be between 0 and 1");
    if (options.disable_z1 && !options.z1_sp_file.empty())
        throw std::runtime_error("--no-z1 cannot be combined with --z1-sp");
    if (options.frame_stride < 1) throw std::runtime_error("frame stride must be positive");
    if (options.time_origin_stride < 1)
        throw std::runtime_error("time-origin stride must be positive");
    if (options.time_origin_count < 2)
        throw std::runtime_error("time-origin count must be at least 2");
    if (options.dw_time_origin_stride < 1)
        throw std::runtime_error(
            "Debye-Waller time-origin stride must be positive");
    if (options.dw_time_origin_count < 2)
        throw std::runtime_error(
            "Debye-Waller time-origin count must be at least 2");
    if (options.time_averaged_msd && options.trajectory_file.empty())
        throw std::runtime_error("--time-averaged-msd requires --trajectory");
    if (options.time_averaged_dw &&
        options.debye_waller_trajectory_file.empty())
        throw std::runtime_error(
            "--time-averaged-dw requires --dw-trajectory");
    if (!(options.diffusion_fit_start_fraction >= 0.0 &&
          options.diffusion_fit_start_fraction < 1.0))
        throw std::runtime_error(
            "diffusion fit start fraction must be at least 0 and less than 1");
    if (std::isfinite(options.debye_waller_time_ps) &&
        !(options.debye_waller_time_ps > 0.0))
        throw std::runtime_error("Debye-Waller time must be positive");
    if (!(options.debye_waller_search_start_ps > 0.0) ||
        !(options.debye_waller_search_end_ps >
          options.debye_waller_search_start_ps))
        throw std::runtime_error(
            "Debye-Waller search times must satisfy 0 < start < end");
    return options;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        const DataFile data = parse_data_file(options.data_file, info);
        const ReducedNetwork network = reduce_network(data, info);
        const int bins = std::max(1, static_cast<int>(
            std::ceil(data.box.lz() / options.bin_width)));
        const std::filesystem::path directory = analysis_directory(
            options.data_file, info, options.output_directory);
        pdms_analysis::create_directory(directory);
        const std::string name = safe_case_name(info.case_name);
        std::string z1_sp_file = options.z1_sp_file;
        if (z1_sp_file.empty() && !options.disable_z1) {
            const std::filesystem::path candidate = directory / "Z1+SP.dat";
            if (std::filesystem::exists(candidate)) z1_sp_file = candidate.string();
        }
        Z1Result z1_result;
        const Z1Result *z1 = nullptr;
        if (!z1_sp_file.empty()) {
            z1_result = read_z1_sp(z1_sp_file, data);
            z1 = &z1_result;
        }
        const std::vector<BinProfile> profile =
            static_profile(data, info, network, bins, z1);
        write_profile(directory / ("network_z_profile." + name + ".tsv"),
                      profile, data, info, z1 != nullptr);
        write_folded_profile(
            directory / ("network_z_profile_folded." + name + ".tsv"),
            profile, data, info, z1 != nullptr);
        write_profile_summary(
            directory / ("network_z_profile_summary." + name + ".tsv"),
            profile, data, info, options, z1 != nullptr);
        const SlidingGrid sliding_grid = make_sliding_grid(data.box, options);
        const std::vector<BinProfile> fine_profile = static_profile(
            data, info, network, sliding_grid.fine_bins, z1);
        const SlidingProfile sliding_profile = make_sliding_profile(
            fine_profile, sliding_grid, info.periodic_z());
        write_profile(
            directory / ("network_z_profile_sliding." + name + ".tsv"),
            sliding_profile.windows, data, info, z1 != nullptr,
            &sliding_profile);
        write_folded_profile(
            directory /
                ("network_z_profile_sliding_folded." + name + ".tsv"),
            sliding_profile.windows, data, info, z1 != nullptr,
            &sliding_profile);
        if (!options.debye_waller_trajectory_file.empty())
            write_debye_waller(
                directory / ("debye_waller_global." + name + ".tsv"),
                directory / ("layer_debye_waller." + name + ".tsv"),
                options, info, data);
        if (options.time_averaged_dw)
            write_time_averaged_debye_waller(
                directory /
                    ("debye_waller_time_averaged." + name + ".tsv"),
                directory /
                    ("layer_debye_waller_time_averaged." + name + ".tsv"),
                options, info, data);
        if (!options.trajectory_file.empty())
            write_layer_dynamics(
                directory / ("layer_dynamics." + name + ".tsv"),
                directory / ("layer_diffusion." + name + ".tsv"),
                options, info, data);
        if (options.time_averaged_msd)
            write_time_averaged_layer_dynamics(
                directory /
                    ("layer_dynamics_time_averaged." + name + ".tsv"),
                directory /
                    ("layer_diffusion_time_averaged." + name + ".tsv"),
                options, info, data);
        write_report(directory / ("profile_report." + name + ".txt"),
                     options, info, data, bins, z1_sp_file);
        std::cout << "Network profiles written to " << directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "network_profile_analyzer: " << error.what() << '\n';
        return 1;
    }
}
