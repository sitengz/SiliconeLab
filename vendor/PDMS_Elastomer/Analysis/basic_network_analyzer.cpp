#include "network_common.hpp"

#include <functional>
#include <iostream>

namespace {
using namespace pdms_analysis;

constexpr double kAvogadroAngstrom = 0.602214076;
constexpr double kBoltzmannKcalPerMolK = 0.00198720425864083;
constexpr double kLammpsRealMvv2e = 2390.05736153349;
constexpr double kBoltzmannSI = 1.380649e-23;

struct Options {
    std::string data_file;
    std::string info_file;
    std::string output_directory;
    int histogram_bins = 40;
    double modulus_temperature_k = -1.0;
};

struct ShapeData {
    Vec3 center;
    double gxx = 0.0, gyy = 0.0, gzz = 0.0;
    double gxy = 0.0, gxz = 0.0, gyz = 0.0;
    std::array<double, 3> eigenvalues{{0.0, 0.0, 0.0}};
    Vec3 principal_axis;
    double rg = 0.0;
    double asphericity = 0.0;
    double acylindricity = 0.0;
    double relative_shape_anisotropy = 0.0;
};

struct StrandProperty {
    const EffectiveStrand *strand = nullptr;
    ShapeData shape;
    double contour_mass = 0.0;
    double ree_xy = 0.0;
    double tortuosity = std::numeric_limits<double>::quiet_NaN();
    Vec3 orientation;
    double p2x = std::numeric_limits<double>::quiet_NaN();
    double p2y = std::numeric_limits<double>::quiet_NaN();
    double p2z = std::numeric_limits<double>::quiet_NaN();
};

struct Statistics {
    long long count = 0;
    double mean = std::numeric_limits<double>::quiet_NaN();
    double standard_deviation = std::numeric_limits<double>::quiet_NaN();
    double rms = std::numeric_limits<double>::quiet_NaN();
    double minimum = std::numeric_limits<double>::quiet_NaN();
    double p05 = std::numeric_limits<double>::quiet_NaN();
    double p25 = std::numeric_limits<double>::quiet_NaN();
    double p50 = std::numeric_limits<double>::quiet_NaN();
    double p75 = std::numeric_limits<double>::quiet_NaN();
    double p95 = std::numeric_limits<double>::quiet_NaN();
    double maximum = std::numeric_limits<double>::quiet_NaN();
};

struct GraphSummary {
    long long active_edges = 0;
    long long dangling = 0;
    long long dangling_loops = 0;
    long long self_loops = 0;
    long long isolated = 0;
    long long components = 0;
    long long largest_component = 0;
    long long largest_component_edges = 0;
    long long largest_component_nodes = 0;
    long long cycle_rank = 0;
    std::map<long long, long long> component_edges;
    std::map<long long, long long> component_nodes;
    std::map<long long, long long> degree_distribution;
};

double percentile(const std::vector<double> &sorted, double fraction) {
    if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double position = fraction * static_cast<double>(sorted.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    return sorted[lower] + (position - static_cast<double>(lower)) *
        (sorted[upper] - sorted[lower]);
}

Statistics statistics(std::vector<double> values) {
    values.erase(std::remove_if(values.begin(), values.end(),
        [](double value) { return !std::isfinite(value); }), values.end());
    Statistics result;
    result.count = static_cast<long long>(values.size());
    if (values.empty()) return result;
    std::sort(values.begin(), values.end());
    double sum = 0.0;
    double sum_squared = 0.0;
    for (double value : values) {
        sum += value;
        sum_squared += value * value;
    }
    result.mean = sum / static_cast<double>(values.size());
    result.rms = std::sqrt(sum_squared / static_cast<double>(values.size()));
    double variance = 0.0;
    for (double value : values) {
        const double delta = value - result.mean;
        variance += delta * delta;
    }
    result.standard_deviation = std::sqrt(
        variance / static_cast<double>(values.size()));
    result.minimum = values.front();
    result.p05 = percentile(values, 0.05);
    result.p25 = percentile(values, 0.25);
    result.p50 = percentile(values, 0.50);
    result.p75 = percentile(values, 0.75);
    result.p95 = percentile(values, 0.95);
    result.maximum = values.back();
    return result;
}

ShapeData shape_data(const std::vector<Vec3> &positions) {
    if (positions.empty()) throw std::runtime_error("cannot calculate Rg of empty path");
    ShapeData shape;
    for (const Vec3 &position : positions) shape.center += position;
    shape.center = (1.0 / static_cast<double>(positions.size())) * shape.center;
    for (const Vec3 &position : positions) {
        const Vec3 delta = position - shape.center;
        shape.gxx += delta.x * delta.x;
        shape.gyy += delta.y * delta.y;
        shape.gzz += delta.z * delta.z;
        shape.gxy += delta.x * delta.y;
        shape.gxz += delta.x * delta.z;
        shape.gyz += delta.y * delta.z;
    }
    const double inverse = 1.0 / static_cast<double>(positions.size());
    shape.gxx *= inverse; shape.gyy *= inverse; shape.gzz *= inverse;
    shape.gxy *= inverse; shape.gxz *= inverse; shape.gyz *= inverse;

    double matrix[3][3] = {
        {shape.gxx, shape.gxy, shape.gxz},
        {shape.gxy, shape.gyy, shape.gyz},
        {shape.gxz, shape.gyz, shape.gzz}
    };
    double vectors[3][3] = {{1.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0},
                            {0.0, 0.0, 1.0}};
    for (int iteration = 0; iteration < 50; ++iteration) {
        int p = 0, q = 1;
        double largest = std::fabs(matrix[0][1]);
        for (const auto& pair : {std::pair<int, int>{0, 2}, {1, 2}}) {
            const double value = std::fabs(matrix[pair.first][pair.second]);
            if (value > largest) { largest = value; p = pair.first; q = pair.second; }
        }
        if (largest < 1.0e-13) break;
        const double app = matrix[p][p];
        const double aqq = matrix[q][q];
        const double apq = matrix[p][q];
        const double angle = 0.5 * std::atan2(2.0 * apq, aqq - app);
        const double cosine = std::cos(angle), sine = std::sin(angle);
        for (int index = 0; index < 3; ++index) {
            if (index == p || index == q) continue;
            const double aip = matrix[index][p];
            const double aiq = matrix[index][q];
            matrix[index][p] = matrix[p][index] =
                cosine * aip - sine * aiq;
            matrix[index][q] = matrix[q][index] =
                sine * aip + cosine * aiq;
        }
        matrix[p][p] = cosine * cosine * app -
            2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q][q] = sine * sine * app +
            2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p][q] = matrix[q][p] = 0.0;
        for (int row = 0; row < 3; ++row) {
            const double vp = vectors[row][p], vq = vectors[row][q];
            vectors[row][p] = cosine * vp - sine * vq;
            vectors[row][q] = sine * vp + cosine * vq;
        }
    }
    std::array<std::pair<double, int>, 3> order{{
        {std::max(0.0, matrix[0][0]), 0},
        {std::max(0.0, matrix[1][1]), 1},
        {std::max(0.0, matrix[2][2]), 2}}};
    std::sort(order.begin(), order.end());
    for (int index = 0; index < 3; ++index)
        shape.eigenvalues[static_cast<std::size_t>(index)] =
            order[static_cast<std::size_t>(index)].first;
    const int principal = order[2].second;
    shape.principal_axis = {
        vectors[0][principal], vectors[1][principal], vectors[2][principal]};
    const double l1 = shape.eigenvalues[0];
    const double l2 = shape.eigenvalues[1];
    const double l3 = shape.eigenvalues[2];
    const double trace = l1 + l2 + l3;
    shape.rg = std::sqrt(std::max(0.0, trace));
    shape.asphericity = l3 - 0.5 * (l1 + l2);
    shape.acylindricity = l2 - l1;
    if (trace > 0.0)
        shape.relative_shape_anisotropy =
            1.0 - 3.0 * (l1 * l2 + l2 * l3 + l3 * l1) / (trace * trace);
    return shape;
}

std::vector<StrandProperty> strand_properties(
    const ReducedNetwork &network, const DataFile &data,
    const ModelInfo &info) {
    std::vector<StrandProperty> result;
    result.reserve(network.strands.size());
    for (const EffectiveStrand &strand : network.strands) {
        StrandProperty property;
        property.strand = &strand;
        const auto positions = unwrapped_path(strand.atoms, data, info);
        property.shape = shape_data(positions);
        for (long long atom : strand.atoms)
            property.contour_mass += data.atoms[static_cast<std::size_t>(atom)].mass;
        property.ree_xy = std::hypot(strand.end_to_end.x, strand.end_to_end.y);
        if (strand.end_to_end_length > 0.0) {
            property.tortuosity = strand.contour_length / strand.end_to_end_length;
            property.orientation =
                (1.0 / strand.end_to_end_length) * strand.end_to_end;
            property.p2x = 0.5 * (3.0 * property.orientation.x * property.orientation.x - 1.0);
            property.p2y = 0.5 * (3.0 * property.orientation.y * property.orientation.y - 1.0);
            property.p2z = 0.5 * (3.0 * property.orientation.z * property.orientation.z - 1.0);
        }
        result.push_back(property);
    }
    return result;
}

GraphSummary graph_summary(const ReducedNetwork &network) {
    GraphSummary summary;
    for (const EffectiveStrand &strand : network.strands) {
        if (strand.status == "active") {
            ++summary.active_edges;
            ++summary.component_edges[strand.graph_component];
        } else if (strand.status == "dangling") ++summary.dangling;
        else if (strand.status == "dangling_loop") ++summary.dangling_loops;
        else if (strand.status == "self_loop") ++summary.self_loops;
        else if (strand.status == "isolated") ++summary.isolated;
    }
    for (std::size_t index = 1; index < network.nodes.size(); ++index) {
        const JunctionNode &node = network.nodes[index];
        ++summary.degree_distribution[node.active_degree];
        if (node.component > 0) ++summary.component_nodes[node.component];
    }
    summary.components = static_cast<long long>(summary.component_nodes.size());
    for (const auto &entry : summary.component_edges) {
        const long long component = entry.first;
        const long long edges = entry.second;
        const long long nodes = summary.component_nodes[component];
        summary.cycle_rank += edges - nodes + 1;
        if (edges > summary.largest_component_edges) {
            summary.largest_component = component;
            summary.largest_component_edges = edges;
            summary.largest_component_nodes = nodes;
        }
    }
    summary.cycle_rank += summary.self_loops;
    return summary;
}

double velocity_temperature(const DataFile &data) {
    if (data.velocities_read != data.declared_atoms || data.declared_atoms < 2)
        return std::numeric_limits<double>::quiet_NaN();
    Vec3 momentum;
    double mass = 0.0;
    double mass_velocity_squared = 0.0;
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        mass += atom.mass;
        momentum += atom.mass * atom.velocity;
        mass_velocity_squared += atom.mass * norm2(atom.velocity);
    }
    const double thermal = mass_velocity_squared - norm2(momentum) / mass;
    const long long dof = 3 * data.declared_atoms - 3;
    return kLammpsRealMvv2e * thermal /
           (static_cast<double>(dof) * kBoltzmannKcalPerMolK);
}

std::array<long long, kComponentCount> component_beads(
    const DataFile &data, const ModelInfo &info) {
    std::array<long long, kComponentCount> counts{{0, 0, 0, 0}};
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        ++counts[static_cast<std::size_t>(
            component_for_molecule(atom.molecule, info))];
    }
    return counts;
}

std::array<double, kComponentCount> component_masses(
    const DataFile &data, const ModelInfo &info) {
    std::array<double, kComponentCount> masses{{0.0, 0.0, 0.0, 0.0}};
    for (long long id = 1; id <= data.declared_atoms; ++id) {
        const Atom &atom = data.atoms[static_cast<std::size_t>(id)];
        masses[static_cast<std::size_t>(
            component_for_molecule(atom.molecule, info))] += atom.mass;
    }
    return masses;
}

std::vector<const StrandProperty *> select_group(
    const std::vector<StrandProperty> &properties,
    const std::string &group, long long largest_component) {
    std::vector<const StrandProperty *> selected;
    for (const StrandProperty &property : properties) {
        const EffectiveStrand &strand = *property.strand;
        bool include = group == "all" || group == strand.status ||
            group == "topology:" + strand.parent_topology;
        if (group == "largest_component")
            include = largest_component > 0 &&
                strand.graph_component == largest_component && strand.status == "active";
        if (include) selected.push_back(&property);
    }
    return selected;
}

using Metric = std::pair<std::string,
    std::function<double(const StrandProperty &)>>;

std::vector<Metric> metrics() {
    return {
        {"contour_beads", [](const StrandProperty &p) { return static_cast<double>(p.strand->atoms.size()); }},
        {"contour_bonds", [](const StrandProperty &p) { return static_cast<double>(p.strand->contour_bonds); }},
        {"contour_mass_g_per_mol", [](const StrandProperty &p) { return p.contour_mass; }},
        {"Lc_A", [](const StrandProperty &p) { return p.strand->contour_length; }},
        {"Ree_A", [](const StrandProperty &p) { return p.strand->end_to_end_length; }},
        {"Ree2_A2", [](const StrandProperty &p) { return norm2(p.strand->end_to_end); }},
        {"Ree_xy_A", [](const StrandProperty &p) { return p.ree_xy; }},
        {"Ree_x2_A2", [](const StrandProperty &p) { return p.strand->end_to_end.x * p.strand->end_to_end.x; }},
        {"Ree_y2_A2", [](const StrandProperty &p) { return p.strand->end_to_end.y * p.strand->end_to_end.y; }},
        {"Ree_z2_A2", [](const StrandProperty &p) { return p.strand->end_to_end.z * p.strand->end_to_end.z; }},
        {"Rg_A", [](const StrandProperty &p) { return p.shape.rg; }},
        {"Rg2_A2", [](const StrandProperty &p) { return p.shape.rg * p.shape.rg; }},
        {"straightness", [](const StrandProperty &p) { return p.strand->straightness; }},
        {"tortuosity", [](const StrandProperty &p) { return p.tortuosity; }},
        {"relative_shape_anisotropy", [](const StrandProperty &p) { return p.shape.relative_shape_anisotropy; }},
        {"P2_x", [](const StrandProperty &p) { return p.p2x; }},
        {"P2_y", [](const StrandProperty &p) { return p.p2y; }},
        {"P2_z", [](const StrandProperty &p) { return p.p2z; }},
        {"Lpp_A", [](const StrandProperty &) { return std::numeric_limits<double>::quiet_NaN(); }}
    };
}

void write_strand_properties(
    const std::filesystem::path &path,
    const std::vector<StrandProperty> &properties) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "strand_id\tparent_molecule\tparent_topology\tstatus\tgraph_component"
        << "\tfirst_node\tsecond_node\tfirst_atom\tsecond_atom\tcontour_beads"
        << "\tcontour_bonds\tcontour_mass_g_per_mol\tLc_A\tRee_x_A\tRee_y_A"
        << "\tRee_z_A\tRee_A\tRee2_A2\tRee_xy_A\tRg_A\tRg2_A2"
        << "\tGxx_A2\tGyy_A2\tGzz_A2\tGxy_A2\tGxz_A2\tGyz_A2"
        << "\tlambda1_A2\tlambda2_A2\tlambda3_A2\tasphericity_A2"
        << "\tacylindricity_A2\trelative_shape_anisotropy\tstraightness"
        << "\ttortuosity\tu_x\tu_y\tu_z\tP2_x\tP2_y\tP2_z"
        << "\tcenter_x_A\tcenter_y_A\tcenter_z_A\tprincipal_x\tprincipal_y"
        << "\tprincipal_z\twinding_x\twinding_y\twinding_z\tLpp_A\tLpp_source\n";
    out << std::setprecision(12);
    for (const StrandProperty &property : properties) {
        const EffectiveStrand &s = *property.strand;
        const ShapeData &g = property.shape;
        out << s.id << '\t' << s.parent_molecule << '\t' << s.parent_topology
            << '\t' << s.status << '\t' << s.graph_component << '\t'
            << s.first_node << '\t' << s.second_node << '\t' << s.first_atom
            << '\t' << s.second_atom << '\t' << s.atoms.size() << '\t'
            << s.contour_bonds << '\t' << property.contour_mass << '\t'
            << s.contour_length << '\t' << s.end_to_end.x << '\t'
            << s.end_to_end.y << '\t' << s.end_to_end.z << '\t'
            << s.end_to_end_length << '\t' << norm2(s.end_to_end) << '\t'
            << property.ree_xy << '\t' << g.rg << '\t' << g.rg * g.rg << '\t'
            << g.gxx << '\t' << g.gyy << '\t' << g.gzz << '\t' << g.gxy
            << '\t' << g.gxz << '\t' << g.gyz << '\t' << g.eigenvalues[0]
            << '\t' << g.eigenvalues[1] << '\t' << g.eigenvalues[2] << '\t'
            << g.asphericity << '\t' << g.acylindricity << '\t'
            << g.relative_shape_anisotropy << '\t' << s.straightness << '\t'
            << property.tortuosity << '\t' << property.orientation.x << '\t'
            << property.orientation.y << '\t' << property.orientation.z << '\t'
            << property.p2x << '\t' << property.p2y << '\t' << property.p2z
            << '\t' << g.center.x << '\t' << g.center.y << '\t' << g.center.z
            << '\t' << g.principal_axis.x << '\t' << g.principal_axis.y << '\t'
            << g.principal_axis.z << '\t' << s.winding[0] << '\t' << s.winding[1]
            << '\t' << s.winding[2] << "\tnan\tnot_available_without_primitive_path_analysis\n";
    }
}

void write_statistics(
    const std::filesystem::path &path,
    const std::vector<StrandProperty> &properties,
    const GraphSummary &graph, const ModelInfo &info) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "group\tmetric\tcount\tmean\tstandard_deviation\trms\tminimum"
        << "\tp05\tp25\tp50\tp75\tp95\tmaximum\n";
    std::vector<std::string> groups = {
        "all", "active", "dangling", "self_loop", "dangling_loop",
        "isolated", "largest_component", "topology:" + info.strand_topology};
    const auto all_metrics = metrics();
    out << std::setprecision(12);
    for (const std::string &group : groups) {
        const auto selected = select_group(properties, group, graph.largest_component);
        for (const Metric &metric : all_metrics) {
            std::vector<double> values;
            values.reserve(selected.size());
            for (const StrandProperty *property : selected)
                values.push_back(metric.second(*property));
            const Statistics stats = statistics(std::move(values));
            out << group << '\t' << metric.first << '\t' << stats.count << '\t'
                << stats.mean << '\t' << stats.standard_deviation << '\t'
                << stats.rms << '\t' << stats.minimum << '\t' << stats.p05
                << '\t' << stats.p25 << '\t' << stats.p50 << '\t' << stats.p75
                << '\t' << stats.p95 << '\t' << stats.maximum << '\n';
        }
    }
}

void write_junctions(
    const std::filesystem::path &path, const ReducedNetwork &network) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "node_id\tkind\trepresentative_molecule\tmember_molecule_count"
        << "\tattachment_count\tactive_degree\tdangling_degree\tcomponent"
        << "\tx_A\ty_A\tz_A\n";
    out << std::setprecision(12);
    for (std::size_t index = 1; index < network.nodes.size(); ++index) {
        const JunctionNode &node = network.nodes[index];
        out << node.id << '\t' << node.kind << '\t'
            << node.representative_molecule << '\t' << node.member_molecules.size()
            << '\t' << node.attachment_atoms.size() << '\t' << node.active_degree
            << '\t' << node.dangling_degree << '\t' << node.component << '\t'
            << node.position.x << '\t' << node.position.y << '\t'
            << node.position.z << '\n';
    }
}

void write_histograms(
    const std::filesystem::path &path,
    const std::vector<StrandProperty> &properties,
    const GraphSummary &graph, int bins) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "group\tmetric\tbin\tlower\tupper\tcenter\tcount\tprobability\n";
    const std::vector<std::string> groups = {"all", "active", "largest_component"};
    const std::vector<Metric> selected_metrics = {
        {"Lc_A", [](const StrandProperty &p) { return p.strand->contour_length; }},
        {"Ree_A", [](const StrandProperty &p) { return p.strand->end_to_end_length; }},
        {"Rg_A", [](const StrandProperty &p) { return p.shape.rg; }},
        {"straightness", [](const StrandProperty &p) { return p.strand->straightness; }},
        {"relative_shape_anisotropy", [](const StrandProperty &p) { return p.shape.relative_shape_anisotropy; }}
    };
    out << std::setprecision(12);
    for (const std::string &group : groups) {
        const auto selected = select_group(properties, group, graph.largest_component);
        for (const Metric &metric : selected_metrics) {
            std::vector<double> values;
            for (const StrandProperty *property : selected) {
                const double value = metric.second(*property);
                if (std::isfinite(value)) values.push_back(value);
            }
            if (values.empty()) continue;
            const auto limits = std::minmax_element(values.begin(), values.end());
            double lower = *limits.first, upper = *limits.second;
            if (upper == lower) upper = lower + std::max(1.0, std::fabs(lower) * 1.0e-6);
            const double width = (upper - lower) / bins;
            std::vector<long long> counts(static_cast<std::size_t>(bins), 0);
            for (double value : values) {
                int bin = static_cast<int>((value - lower) / width);
                if (bin == bins) --bin;
                ++counts[static_cast<std::size_t>(bin)];
            }
            for (int bin = 0; bin < bins; ++bin) {
                const double bin_lower = lower + bin * width;
                out << group << '\t' << metric.first << '\t' << bin + 1 << '\t'
                    << bin_lower << '\t' << bin_lower + width << '\t'
                    << bin_lower + 0.5 * width << '\t'
                    << counts[static_cast<std::size_t>(bin)] << '\t'
                    << static_cast<double>(counts[static_cast<std::size_t>(bin)]) /
                       static_cast<double>(values.size()) << '\n';
            }
        }
    }
}

void write_network_statistics(
    const std::filesystem::path &path, const DataFile &data,
    const ModelInfo &info,
    const ReducedNetwork &network,
    const GraphSummary &graph, const std::vector<StrandProperty> &properties,
    double modulus_temperature) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "property\tvalue\tunits\tdefinition\n" << std::setprecision(12);
    const double volume = data.box.lx() * data.box.ly() * data.box.lz();
    const double density = data.total_mass_g_per_mol /
        (kAvogadroAngstrom * volume);
    const double active_density = graph.active_edges / volume;
    const double affine_modulus_mpa =
        active_density * 1.0e30 * kBoltzmannSI * modulus_temperature / 1.0e6;
    double degree_sum = 0.0;
    long long degree_nodes = 0;
    for (std::size_t index = 1; index < network.nodes.size(); ++index) {
        const long long degree = network.nodes[index].active_degree;
        if (degree > 0) { degree_sum += degree; ++degree_nodes; }
    }
    const double mean_functionality = degree_nodes > 0
        ? degree_sum / degree_nodes : std::numeric_limits<double>::quiet_NaN();
    const double phantom_factor = std::isfinite(mean_functionality) && mean_functionality > 2.0
        ? 1.0 - 2.0 / mean_functionality : std::numeric_limits<double>::quiet_NaN();
    const auto active = select_group(properties, "active", graph.largest_component);
    std::vector<double> masses;
    for (const StrandProperty *property : active) masses.push_back(property->contour_mass);
    const Statistics mass_stats = statistics(std::move(masses));
    const auto beads = component_beads(data, info);
    const auto component_mass = component_masses(data, info);
    long long crosslink_bonds = 0;
    for (const Bond &bond : data.bonds)
        if (bond.type == info.crosslink_bond_type) ++crosslink_bonds;
    const auto line = [&](const std::string &key, double value,
                          const std::string &units, const std::string &definition) {
        out << key << '\t' << value << '\t' << units << '\t' << definition << '\n';
    };
    line("atoms", static_cast<double>(data.declared_atoms), "count", "LAMMPS header");
    line("molecules", static_cast<double>(info.total_molecules), "count", "info-file component total");
    line("bonds", static_cast<double>(data.declared_bonds), "count", "LAMMPS header");
    line("angles", static_cast<double>(data.declared_angles), "count", "LAMMPS header");
    line("dihedrals", static_cast<double>(data.declared_dihedrals), "count", "LAMMPS header");
    line("volume", volume, "A^3", "orthorhombic box volume");
    line("box_Lx", data.box.lx(), "A", "complete x box length");
    line("box_Ly", data.box.ly(), "A", "complete y box length");
    line("box_Lz", data.box.lz(), "A", "complete z box length");
    line("total_mass", data.total_mass_g_per_mol, "g/mol-equivalent", "sum of bead masses");
    line("density", density, "g/cm^3", "mass divided by box volume");
    line("velocities_read", static_cast<double>(data.velocities_read), "count", "parsed velocity records");
    line("velocity_temperature", velocity_temperature(data), "K", "COM-corrected kinetic temperature; NaN without complete velocities");
    const std::array<std::string, kComponentCount> component_names{{
        "strands", "crosslinkers", "moderators", "filler"}};
    for (int component = 0; component < kComponentCount; ++component) {
        const std::size_t index = static_cast<std::size_t>(component);
        const std::string prefix = "component_" + component_names[index] + "_";
        line(prefix + "molecules",
             static_cast<double>(info.components[index].molecules), "count",
             "info-file molecule count");
        line(prefix + "beads", static_cast<double>(beads[index]), "count",
             "snapshot bead count");
        line(prefix + "mass", component_mass[index], "g/mol-equivalent",
             "snapshot bead-mass sum");
        line(prefix + "mass_fraction", component_mass[index] /
             data.total_mass_g_per_mol, "fraction", "component mass / total mass");
    }
    line("crosslink_bonds", static_cast<double>(crosslink_bonds), "count",
         "bonds with the info-file crosslink bond type");
    line("strand_sites_total", static_cast<double>(network.strand_reactive_total),
         "count", "component-1 functional sites");
    line("strand_sites_reacted", static_cast<double>(network.strand_reactive_reacted),
         "count", "reacted component-1 functional sites");
    line("strand_site_conversion", network.strand_reactive_total == 0 ? 0.0 :
         static_cast<double>(network.strand_reactive_reacted) / network.strand_reactive_total,
         "fraction", "reacted component-1 functional sites");
    line("crosslinker_sites_total",
         static_cast<double>(network.crosslinker_reactive_total), "count",
         "component-2 functional sites");
    line("crosslinker_sites_reacted",
         static_cast<double>(network.crosslinker_reactive_reacted), "count",
         "reacted component-2 functional sites");
    line("crosslinker_site_conversion", network.crosslinker_reactive_total == 0 ? 0.0 :
         static_cast<double>(network.crosslinker_reactive_reacted) / network.crosslinker_reactive_total,
         "fraction", "reacted component-2 functional sites");
    line("junction_nodes", static_cast<double>(network.nodes.size() - 1), "count", "reduced graph nodes");
    line("effective_strands", static_cast<double>(network.strands.size()), "count", "all reduced paths");
    line("active_strands", static_cast<double>(graph.active_edges), "count", "edges between distinct junctions");
    line("dangling_strands", static_cast<double>(graph.dangling), "count", "one-junction effective strands");
    line("dangling_loops", static_cast<double>(graph.dangling_loops), "count", "rings with exactly one reacted site");
    line("self_loops", static_cast<double>(graph.self_loops), "count", "effective strands whose endpoints share a junction");
    line("isolated_strands", static_cast<double>(graph.isolated), "count", "effective strands with no junction endpoint");
    line("active_strand_density", active_density, "A^-3", "active edges divided by box volume");
    line("active_graph_components", static_cast<double>(graph.components), "count", "connected components containing active junctions");
    line("largest_component_nodes", static_cast<double>(graph.largest_component_nodes), "count", "junctions in largest active component");
    line("largest_component_edges", static_cast<double>(graph.largest_component_edges), "count", "active strands in largest component");
    line("largest_component_active_fraction", graph.active_edges == 0 ? 0.0 :
         static_cast<double>(graph.largest_component_edges) / graph.active_edges,
         "fraction", "active edges in largest graph component");
    line("cycle_rank", static_cast<double>(graph.cycle_rank), "count", "sum(E-V+1) plus self-loops");
    line("mean_active_junction_functionality", mean_functionality, "degree", "mean degree over nonzero-degree nodes");
    line("mean_active_contour_mass", mass_stats.mean, "g/mol", "contour-mass proxy; shared endpoints may repeat");
    line("affine_modulus_estimate", affine_modulus_mpa, "MPa", "nu_active k_B T structural estimate");
    line("phantom_modulus_estimate", affine_modulus_mpa * phantom_factor, "MPa", "affine estimate times (1-2/f_mean)");
    line("modulus_temperature", modulus_temperature, "K", "temperature used for structural modulus estimates");
    line("Lpp_available", 0.0, "boolean", "primitive-path results deliberately not inferred from snapshot");
}

void write_report(
    const std::filesystem::path &path, const Options &options,
    const DataFile &data, const ModelInfo &info,
    const ReducedNetwork &network, const GraphSummary &graph,
    const std::vector<StrandProperty> &properties,
    double modulus_temperature) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    const double volume = data.box.lx() * data.box.ly() * data.box.lz();
    const auto beads = component_beads(data, info);
    const auto masses = component_masses(data, info);
    const double temperature = velocity_temperature(data);
    std::map<std::string, long long> statuses;
    for (const EffectiveStrand &strand : network.strands) ++statuses[strand.status];
    const auto active = select_group(properties, "active", graph.largest_component);
    std::vector<double> ree, rg, lc, straightness;
    for (const StrandProperty *property : active) {
        ree.push_back(property->strand->end_to_end_length);
        rg.push_back(property->shape.rg);
        lc.push_back(property->strand->contour_length);
        straightness.push_back(property->strand->straightness);
    }
    const Statistics ree_stats = statistics(ree);
    const Statistics rg_stats = statistics(rg);
    const Statistics lc_stats = statistics(lc);
    const Statistics straight_stats = statistics(straightness);
    long long crosslink_bonds = 0;
    for (const Bond &bond : data.bonds)
        if (bond.type == info.crosslink_bond_type) ++crosslink_bonds;

    out << std::setprecision(10)
        << "PDMS basic network information report\n"
        << "case: " << info.case_name << "\n"
        << "source: " << options.data_file << "\n"
        << "geometry: " << info.geometry << "\n"
        << "strand topology: " << info.strand_topology << "\n\n"
        << "System\n"
        << "  atoms / bonds / angles / dihedrals: " << data.declared_atoms << " / "
        << data.declared_bonds << " / " << data.declared_angles << " / "
        << data.declared_dihedrals << "\n"
        << "  box Lx Ly Lz (A): " << data.box.lx() << ' ' << data.box.ly() << ' '
        << data.box.lz() << "\n"
        << "  volume (A^3): " << volume << "\n"
        << "  total mass (g/mol-equivalent): " << data.total_mass_g_per_mol << "\n"
        << "  density (g/cm^3): " << data.total_mass_g_per_mol /
            (kAvogadroAngstrom * volume) << "\n"
        << "  velocity temperature (K): " << temperature << "\n"
        << "  target final temperature (K): " << info.final_temperature_k << "\n";
    const std::array<std::string, kComponentCount> names{{
        "strands", "crosslinkers", "moderators", "filler"}};
    out << "\nComponents\n";
    for (int component = 0; component < kComponentCount; ++component) {
        const std::size_t index = static_cast<std::size_t>(component);
        out << "  " << names[index] << ": molecules "
            << info.components[index].molecules << ", beads " << beads[index]
            << ", mass " << masses[index] << " g/mol-equivalent, mass fraction "
            << (data.total_mass_g_per_mol == 0.0 ? 0.0 :
                masses[index] / data.total_mass_g_per_mol) << "\n";
    }
    out << "\nReaction and graph\n"
        << "  crosslink bonds: " << crosslink_bonds << "\n"
        << "  strand-site conversion: " << network.strand_reactive_reacted << " / "
        << network.strand_reactive_total << "\n"
        << "  crosslinker-site conversion: " << network.crosslinker_reactive_reacted
        << " / " << network.crosslinker_reactive_total << "\n"
        << "  junction nodes: " << network.nodes.size() - 1 << "\n"
        << "  effective strands: " << network.strands.size() << "\n"
        << "  active components: " << graph.components << "\n"
        << "  largest component edges / nodes: " << graph.largest_component_edges
        << " / " << graph.largest_component_nodes << "\n"
        << "  cycle rank: " << graph.cycle_rank << "\n";
    for (const auto &entry : statuses)
        out << "  " << entry.first << ": " << entry.second << "\n";
    out << "  degree distribution (degree:nodes):";
    for (const auto &entry : graph.degree_distribution)
        out << ' ' << entry.first << ':' << entry.second;
    out << "\n\nActive-strand conformation\n"
        << "  count: " << active.size() << "\n"
        << "  mean / RMS Ree (A): " << ree_stats.mean << " / " << ree_stats.rms << "\n"
        << "  Ree p05/p25/p50/p75/p95 (A): " << ree_stats.p05 << ' '
        << ree_stats.p25 << ' ' << ree_stats.p50 << ' ' << ree_stats.p75 << ' '
        << ree_stats.p95 << "\n"
        << "  mean / RMS Rg (A): " << rg_stats.mean << " / " << rg_stats.rms << "\n"
        << "  mean Lc (A): " << lc_stats.mean << "\n"
        << "  mean straightness Ree/Lc: " << straight_stats.mean << "\n"
        << "\nStructural elasticity estimates\n"
        << "  temperature (K): " << modulus_temperature << "\n"
        << "  estimates are reported in network_statistics and are not measured moduli\n"
        << "\nPrimitive path\n"
        << "  Lpp: unavailable (NaN)\n"
        << "  Reason: Lpp requires validated primitive-path coordinates/results; "
           "Lc, Ree, and graph shortest paths are not substitutes.\n";
}

void print_help(const char *program) {
    std::cout
        << "Usage: " << program << " <case>.npt_eq <case>.info [options]\n\n"
        << "Options:\n"
        << "  --output-dir PATH          output directory (default analysis_<case>)\n"
        << "  --histogram-bins N         number of bins (default 40)\n"
        << "  --modulus-temperature X    structural estimate temperature in K\n"
        << "  --help                     show this help\n";
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
        if (option == "--output-dir") options.output_directory = value();
        else if (option == "--histogram-bins") options.histogram_bins = std::stoi(value());
        else if (option == "--modulus-temperature")
            options.modulus_temperature_k = std::stod(value());
        else if (option == "--help") { print_help(argv[0]); std::exit(0); }
        else throw std::runtime_error("unknown option: " + option);
    }
    if (options.histogram_bins < 1 || options.histogram_bins > 10000)
        throw std::runtime_error("histogram bins must be between 1 and 10000");
    if (options.modulus_temperature_k == 0.0 ||
        options.modulus_temperature_k < -1.0)
        throw std::runtime_error("modulus temperature must be positive");
    return options;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        const DataFile data = parse_data_file(options.data_file, info);
        if (data.declared_atoms != info.total_beads)
            throw std::runtime_error(
                "snapshot atom count does not match info-file total_beads");
        for (const ComponentInfo &component : info.components) {
            for (long long molecule = component.molecule_start;
                 molecule <= component.molecule_end && component.molecules > 0;
                 ++molecule) {
                const auto &atoms = data.molecule_atoms[
                    static_cast<std::size_t>(molecule)];
                if (static_cast<long long>(atoms.size()) !=
                    component.beads_per_molecule)
                    throw std::runtime_error(
                        "snapshot molecule " + std::to_string(molecule) +
                        " bead count does not match info file");
            }
        }
        const ReducedNetwork network = reduce_network(data, info);
        const std::vector<StrandProperty> properties =
            strand_properties(network, data, info);
        const GraphSummary graph = graph_summary(network);
        const double modulus_temperature = options.modulus_temperature_k > 0.0
            ? options.modulus_temperature_k : info.final_temperature_k;
        const std::filesystem::path directory = analysis_directory(
            options.data_file, info, options.output_directory);
        pdms_analysis::create_directory(directory);
        const std::string name = safe_case_name(info.case_name);
        write_strand_properties(
            directory / ("strand_properties." + name + ".tsv"), properties);
        write_statistics(directory / ("strand_statistics." + name + ".tsv"),
                         properties, graph, info);
        write_junctions(directory / ("junction_properties." + name + ".tsv"),
                        network);
        write_histograms(directory / ("strand_histograms." + name + ".tsv"),
                         properties, graph, options.histogram_bins);
        write_network_statistics(
            directory / ("network_statistics." + name + ".tsv"),
            data, info, network, graph, properties, modulus_temperature);
        write_report(directory / ("basic_network_report." + name + ".txt"),
                     options, data, info, network, graph, properties,
                     modulus_temperature);
        std::cout << "Basic network information written to "
                  << directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "basic_network_analyzer: " << error.what() << '\n';
        return 1;
    }
}
