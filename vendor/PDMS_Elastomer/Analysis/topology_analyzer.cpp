#include "network_common.hpp"

#include <functional>
#include <iostream>

namespace {
using namespace pdms_analysis;

struct Options {
    std::string data_file;
    std::string info_file;
    std::string output_directory;
    std::string z1_selection = "network";
    double z1_max_bond = 0.0;
    bool z1_scaling_requested = false;
    int image_search_bound = 2;
    bool skip_self_paths = false;
};

struct DirectionalPath {
    long long source_node = 0;
    int axis = 0;
    double contour_length = std::numeric_limits<double>::infinity();
    long long strand_count = -1;
};

struct GraphArc {
    long long neighbor = 0;
    std::array<int, 3> shift{{0, 0, 0}};
    double weight = 0.0;
};

struct Z1Chain {
    long long effective_strand_id = 0;
    long long parent_molecule = 0;
    std::string parent_topology;
    std::string status;
    long long first_node = 0;
    long long second_node = 0;
    long long first_atom = 0;
    long long second_atom = 0;
    long long effective_contour_beads = 0;
    long long reacted_site_index = 0;
    long long next_graft_atom = 0;
    std::vector<long long> atoms;
    std::vector<long long> excluded_atoms;
};

struct Z1Summary {
    double coordinate_scale = 1.0;
    long long selected_paths = 0;
    long long chains_written = 0;
    long long empty_strands_skipped = 0;
    long long excluded_contour_occurrences = 0;
};

std::vector<std::vector<GraphArc>> graph_arcs(
    const ReducedNetwork &network) {
    std::vector<std::vector<GraphArc>> adjacency(network.nodes.size());
    for (const EffectiveStrand &strand : network.strands) {
        if (strand.status != "active") continue;
        adjacency[static_cast<std::size_t>(strand.first_node)].push_back({
            strand.second_node, strand.winding, strand.contour_length});
        std::array<int, 3> reverse{{
            -strand.winding[0], -strand.winding[1], -strand.winding[2]}};
        adjacency[static_cast<std::size_t>(strand.second_node)].push_back({
            strand.first_node, reverse, strand.contour_length});
    }
    return adjacency;
}

long long shift_index(
    const std::array<int, 3> &shift, int dimensions, int bound) {
    const int side = 2 * bound + 1;
    long long index = 0;
    long long multiplier = 1;
    for (int axis = 0; axis < dimensions; ++axis) {
        if (shift[axis] < -bound || shift[axis] > bound) return -1;
        index += multiplier * (shift[axis] + bound);
        multiplier *= side;
    }
    return index;
}

std::vector<DirectionalPath> directional_self_paths(
    const ReducedNetwork &network, bool periodic_z, int bound) {
    const int dimensions = periodic_z ? 3 : 2;
    const int side = 2 * bound + 1;
    long long shift_states = 1;
    for (int axis = 0; axis < dimensions; ++axis) shift_states *= side;
    const auto adjacency = graph_arcs(network);
    std::vector<DirectionalPath> results;
    for (std::size_t source = 1; source < network.nodes.size(); ++source) {
        if (adjacency[source].empty()) continue;
        const std::size_t state_count = network.nodes.size() *
            static_cast<std::size_t>(shift_states);
        std::vector<double> distance(state_count,
            std::numeric_limits<double>::infinity());
        std::vector<long long> hops(state_count,
            std::numeric_limits<long long>::max());
        const std::array<int, 3> zero{{0, 0, 0}};
        const long long zero_shift = shift_index(zero, dimensions, bound);
        const std::size_t start = source * static_cast<std::size_t>(shift_states) +
            static_cast<std::size_t>(zero_shift);
        using QueueEntry = std::tuple<double, long long, long long,
                                      std::array<int, 3>>;
        std::priority_queue<QueueEntry, std::vector<QueueEntry>,
                            std::greater<QueueEntry>> queue;
        distance[start] = 0.0;
        hops[start] = 0;
        queue.push({0.0, 0, static_cast<long long>(source), zero});
        while (!queue.empty()) {
            const auto [current_distance, current_hops, node, shift] = queue.top();
            queue.pop();
            const long long current_shift = shift_index(shift, dimensions, bound);
            const std::size_t current_state =
                static_cast<std::size_t>(node) * static_cast<std::size_t>(shift_states) +
                static_cast<std::size_t>(current_shift);
            if (current_distance > distance[current_state] + 1.0e-12 ||
                current_hops != hops[current_state]) continue;
            for (const GraphArc &arc : adjacency[static_cast<std::size_t>(node)]) {
                std::array<int, 3> next_shift = shift;
                for (int axis = 0; axis < dimensions; ++axis)
                    next_shift[axis] += arc.shift[axis];
                const long long next_shift_id =
                    shift_index(next_shift, dimensions, bound);
                if (next_shift_id < 0) continue;
                const std::size_t next_state =
                    static_cast<std::size_t>(arc.neighbor) *
                    static_cast<std::size_t>(shift_states) +
                    static_cast<std::size_t>(next_shift_id);
                const double candidate = current_distance + arc.weight;
                const long long candidate_hops = current_hops + 1;
                if (candidate + 1.0e-12 < distance[next_state] ||
                    (std::fabs(candidate - distance[next_state]) <= 1.0e-12 &&
                     candidate_hops < hops[next_state])) {
                    distance[next_state] = candidate;
                    hops[next_state] = candidate_hops;
                    queue.push({candidate, candidate_hops, arc.neighbor, next_shift});
                }
            }
        }
        for (int axis = 0; axis < dimensions; ++axis) {
            std::array<int, 3> target{{0, 0, 0}};
            target[axis] = 1;
            const long long target_shift = shift_index(target, dimensions, bound);
            const std::size_t target_state =
                source * static_cast<std::size_t>(shift_states) +
                static_cast<std::size_t>(target_shift);
            results.push_back({static_cast<long long>(source), axis,
                               distance[target_state], hops[target_state] ==
                                   std::numeric_limits<long long>::max()
                                   ? -1 : hops[target_state]});
        }
    }
    return results;
}

double percentile(std::vector<double> values, double fraction) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double position = fraction * (values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    return values[lower] + (position - lower) * (values[upper] - values[lower]);
}

std::string axis_name(int axis) {
    return axis == 0 ? "x" : axis == 1 ? "y" : "z";
}

bool select_for_z1(
    const std::string &status, bool parent_reacted,
    const std::string &selection) {
    if (selection == "network")
        return parent_reacted && status != "isolated";
    if (selection == "active") return status == "active";
    if (selection == "active-and-dangling")
        return parent_reacted && (status == "active" ||
            status == "dangling" || status == "dangling_loop");
    if (selection == "all-linearizable") return true;
    throw std::runtime_error("unknown Z1 selection: " + selection);
}

void append_unique(std::vector<long long> &values, long long value) {
    if (std::find(values.begin(), values.end(), value) == values.end())
        values.push_back(value);
}

Z1Chain z1_chain_metadata(const EffectiveStrand &strand) {
    Z1Chain chain;
    chain.effective_strand_id = strand.id;
    chain.parent_molecule = strand.parent_molecule;
    chain.parent_topology = strand.parent_topology;
    chain.status = strand.status;
    chain.first_node = strand.first_node;
    chain.second_node = strand.second_node;
    chain.first_atom = strand.first_atom;
    chain.second_atom = strand.second_atom;
    chain.effective_contour_beads =
        static_cast<long long>(strand.atoms.size());
    return chain;
}

Z1Chain z1_chain(
    const EffectiveStrand &strand, const DataFile &data,
    const ModelInfo &info) {
    Z1Chain chain = z1_chain_metadata(strand);
    if (strand.parent_topology == "ring") {
        if (!strand.atoms.empty()) {
            append_unique(chain.excluded_atoms, strand.atoms.front());
            append_unique(chain.excluded_atoms, strand.atoms.back());
        }
        if (strand.atoms.size() > 2)
            chain.atoms.assign(strand.atoms.begin() + 1, strand.atoms.end() - 1);
    } else if (strand.parent_topology == "star") {
        const auto &molecule_atoms = data.molecule_atoms[
            static_cast<std::size_t>(strand.parent_molecule)];
        for (long long atom : strand.atoms) {
            if (local_rank(atom, molecule_atoms) <= info.star_center_count)
                append_unique(chain.excluded_atoms, atom);
            else
                chain.atoms.push_back(atom);
        }
    } else {
        // Linear paths retain their complete reduced contours.
        chain.atoms = strand.atoms;
    }
    return chain;
}

long long graft_atom_for_site(
    long long site, const std::vector<long long> &molecule_atoms,
    const ModelInfo &info) {
    const long long rank = local_rank(site, molecule_atoms);
    const long long offset = rank - info.graft_backbone_length;
    if (info.graft_backbone_length < 1 || info.graft_side_chain_length < 1 ||
        offset < 1 || offset % info.graft_side_chain_length != 0)
        throw std::runtime_error(
            "reacted grafted site is not a generated side-chain end");
    const long long side = offset / info.graft_side_chain_length - 1;
    if (side < 0 || side >= info.graft_side_chain_count)
        throw std::runtime_error("reacted grafted site has invalid side-chain index");
    const long long graft_rank = 1 + side * (info.graft_spacing + 1);
    if (graft_rank < 1 || graft_rank > info.graft_backbone_length)
        throw std::runtime_error("grafted side chain has invalid backbone attachment");
    return molecule_atoms[static_cast<std::size_t>(graft_rank - 1)];
}

std::vector<Z1Chain> grafted_z1_chains(
    const ReducedNetwork &network, const DataFile &data,
    const ModelInfo &info, const Options &options) {
    const auto internal = strand_internal_adjacency(data, info);
    std::map<std::pair<long long, long long>, const EffectiveStrand *> segments;
    for (const EffectiveStrand &strand : network.strands) {
        if (strand.parent_topology != "grafted") continue;
        segments[{strand.first_atom, strand.second_atom}] = &strand;
    }

    std::vector<Z1Chain> result;
    const ComponentInfo &component = info.components[kStrand];
    for (long long molecule = component.molecule_start;
         molecule > 0 && molecule <= component.molecule_end; ++molecule) {
        const auto &molecule_atoms = data.molecule_atoms[
            static_cast<std::size_t>(molecule)];
        std::vector<long long> reacted_sites;
        for (const auto &entry : network.network_atom_node) {
            const long long atom = entry.first;
            if (data.atoms[static_cast<std::size_t>(atom)].molecule == molecule)
                reacted_sites.push_back(atom);
        }
        std::sort(reacted_sites.begin(), reacted_sites.end(),
            [&](long long first, long long second) {
                return graft_atom_for_site(first, molecule_atoms, info) <
                    graft_atom_for_site(second, molecule_atoms, info);
            });
        if (reacted_sites.empty()) continue;

        for (std::size_t site = 0; site < reacted_sites.size(); ++site) {
            const long long current = reacted_sites[site];
            const long long current_graft =
                graft_atom_for_site(current, molecule_atoms, info);
            Z1Chain chain;
            if (site + 1 < reacted_sites.size()) {
                const long long next = reacted_sites[site + 1];
                const long long next_graft =
                    graft_atom_for_site(next, molecule_atoms, info);
                const auto found = segments.find({current, next});
                if (found == segments.end())
                    throw std::runtime_error(
                        "missing reduced grafted segment between ordered reacted sites");
                chain = z1_chain_metadata(*found->second);
                chain.reacted_site_index = static_cast<long long>(site + 1);
                chain.next_graft_atom = next_graft;
                const std::vector<long long> path =
                    unique_path(current, next, internal);
                const auto stop = std::find(path.begin(), path.end(), next_graft);
                if (stop == path.end() || stop == path.begin())
                    throw std::runtime_error(
                        "next grafting point is not internal to grafted path");
                chain.atoms.assign(path.begin(), stop);
                chain.excluded_atoms.assign(stop, path.end());
                chain.first_atom = chain.atoms.front();
                chain.second_atom = chain.atoms.back();
            } else {
                std::vector<long long> path =
                    unique_path(current, current_graft, internal);
                chain.parent_molecule = molecule;
                chain.parent_topology = "grafted";
                chain.status = "dangling";
                chain.first_node = network.network_atom_node.at(current);
                chain.first_atom = current;
                chain.effective_contour_beads =
                    static_cast<long long>(path.size());
                chain.reacted_site_index = static_cast<long long>(site + 1);
                chain.excluded_atoms.push_back(current_graft);
                path.pop_back();
                if (path.empty())
                    throw std::runtime_error(
                        "final grafted side chain has no side-chain beads");
                chain.atoms = path;
                chain.second_atom = chain.atoms.back();
            }
            if (select_for_z1(chain.status, true, options.z1_selection))
                result.push_back(std::move(chain));
        }
    }
    return result;
}

void write_strands(
    const std::filesystem::path &path, const ReducedNetwork &network) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "strand_id\tparent_molecule\tparent_topology\tstatus\tfirst_node"
        << "\tsecond_node\tfirst_atom\tsecond_atom\tcontour_beads"
        << "\tcontour_bonds\tcontour_length_A\tRee_x_A\tRee_y_A\tRee_z_A"
        << "\tRee_A\tRee2_A2\tstraightness\twinding_x\twinding_y"
        << "\twinding_z\tgraph_component\n";
    out << std::setprecision(12);
    for (const EffectiveStrand &strand : network.strands) {
        out << strand.id << '\t' << strand.parent_molecule << '\t'
            << strand.parent_topology << '\t' << strand.status << '\t'
            << strand.first_node << '\t' << strand.second_node << '\t'
            << strand.first_atom << '\t' << strand.second_atom << '\t'
            << strand.atoms.size() << '\t' << strand.contour_bonds << '\t'
            << strand.contour_length << '\t' << strand.end_to_end.x << '\t'
            << strand.end_to_end.y << '\t' << strand.end_to_end.z << '\t'
            << strand.end_to_end_length << '\t'
            << strand.end_to_end_length * strand.end_to_end_length << '\t'
            << strand.straightness << '\t' << strand.winding[0] << '\t'
            << strand.winding[1] << '\t' << strand.winding[2] << '\t'
            << strand.graph_component << '\n';
    }
}

void write_nodes(
    const std::filesystem::path &path, const ReducedNetwork &network) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "node_id\tkind\trepresentative_molecule\tmember_molecules"
        << "\tattachments\tactive_degree\tdangling_degree\tcomponent"
        << "\tx_A\ty_A\tz_A\n";
    out << std::setprecision(12);
    for (std::size_t index = 1; index < network.nodes.size(); ++index) {
        const JunctionNode &node = network.nodes[index];
        out << node.id << '\t' << node.kind << '\t'
            << node.representative_molecule << '\t';
        for (std::size_t i = 0; i < node.member_molecules.size(); ++i) {
            if (i) out << ',';
            out << node.member_molecules[i];
        }
        out << '\t' << node.attachment_atoms.size() << '\t'
            << node.active_degree << '\t' << node.dangling_degree << '\t'
            << node.component << '\t' << node.position.x << '\t'
            << node.position.y << '\t' << node.position.z << '\n';
    }
}

void write_parents(
    const std::filesystem::path &path, const ReducedNetwork &network) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "molecule\ttopology\tpossible_sites\treacted_sites"
        << "\teffective_strands\tstate\n";
    for (const ParentRecord &parent : network.parents)
        out << parent.molecule << '\t' << parent.topology << '\t'
            << parent.possible_sites << '\t' << parent.reacted_sites << '\t'
            << parent.effective_strands << '\t' << parent.state << '\n';
}

void write_directional_paths(
    const std::filesystem::path &path,
    const std::vector<DirectionalPath> &paths) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << "source_node\tdirection\treachable\tshortest_contour_length_A"
        << "\tstrand_count\n";
    out << std::setprecision(12);
    for (const DirectionalPath &entry : paths) {
        out << entry.source_node << '\t' << axis_name(entry.axis) << '\t'
            << (std::isfinite(entry.contour_length) ? 1 : 0) << '\t';
        if (std::isfinite(entry.contour_length)) out << entry.contour_length;
        else out << "nan";
        out << '\t' << entry.strand_count << '\n';
    }
}

double z1_scale(
    const std::vector<Z1Chain> &selected,
    const DataFile &data, const ModelInfo &info, double target) {
    double maximum = 0.0;
    for (const Z1Chain &chain : selected) {
        const auto positions = unwrapped_path(chain.atoms, data, info);
        for (std::size_t i = 1; i < positions.size(); ++i)
            maximum = std::max(maximum, norm(positions[i] - positions[i - 1]));
    }
    return maximum > target ? target / maximum : 1.0;
}

Z1Summary write_z1(
    const std::filesystem::path &config_path,
    const std::filesystem::path &map_path,
    const ReducedNetwork &network, const DataFile &data,
    const ModelInfo &info, const Options &options) {
    Z1Summary summary;
    std::set<long long> reacted_parents;
    for (const ParentRecord &parent : network.parents)
        if (parent.reacted_sites > 0) reacted_parents.insert(parent.molecule);
    std::vector<Z1Chain> selected;
    if (info.strand_topology == "grafted") {
        selected = grafted_z1_chains(network, data, info, options);
        summary.selected_paths = static_cast<long long>(selected.size());
        for (const Z1Chain &chain : selected)
            summary.excluded_contour_occurrences +=
                static_cast<long long>(chain.excluded_atoms.size());
    } else {
        for (const EffectiveStrand &strand : network.strands) {
            const bool parent_reacted =
                reacted_parents.count(strand.parent_molecule) != 0;
            if (!select_for_z1(
                    strand.status, parent_reacted, options.z1_selection))
                continue;
            ++summary.selected_paths;
            Z1Chain chain = z1_chain(strand, data, info);
            summary.excluded_contour_occurrences +=
                static_cast<long long>(strand.atoms.size() - chain.atoms.size());
            if (chain.atoms.empty()) {
                ++summary.empty_strands_skipped;
                continue;
            }
            selected.push_back(std::move(chain));
        }
    }
    summary.chains_written = static_cast<long long>(selected.size());

    std::map<long long, long long> atom_owner;
    for (std::size_t chain = 0; chain < selected.size(); ++chain) {
        for (long long atom : selected[chain].atoms) {
            const auto inserted = atom_owner.emplace(
                atom, static_cast<long long>(chain + 1));
            if (!inserted.second)
                throw std::runtime_error(
                    "Z1 bead overlap remains between chains " +
                    std::to_string(inserted.first->second) + " and " +
                    std::to_string(chain + 1) + " at atom " +
                    std::to_string(atom));
        }
    }

    summary.coordinate_scale = options.z1_scaling_requested
        ? z1_scale(selected, data, info, options.z1_max_bond) : 1.0;
    std::ofstream config(config_path);
    if (!config) throw std::runtime_error("cannot write " + config_path.string());
    config << selected.size() << '\n' << std::setprecision(15)
           << summary.coordinate_scale * data.box.lx() << ' '
           << summary.coordinate_scale * data.box.ly() << ' '
           << summary.coordinate_scale * data.box.lz() << '\n';
    for (std::size_t i = 0; i < selected.size(); ++i) {
        if (i) config << ' ';
        config << selected[i].atoms.size();
    }
    config << '\n';
    for (const Z1Chain &chain : selected) {
        for (const Vec3 &position : unwrapped_path(chain.atoms, data, info))
            config << summary.coordinate_scale * position.x << ' '
                   << summary.coordinate_scale * position.y << ' '
                   << summary.coordinate_scale * position.z << '\n';
    }

    std::ofstream mapping(map_path);
    if (!mapping) throw std::runtime_error("cannot write " + map_path.string());
    mapping << "z1_chain_id\teffective_strand_id\tparent_molecule"
            << "\tparent_topology\tstatus\tfirst_node\tsecond_node"
            << "\tfirst_atom\tsecond_atom\teffective_contour_beads\tz1_beads"
            << "\treacted_site_index\tnext_graft_atom\texcluded_atom_ids"
            << "\tz1_atom_ids\tcoordinate_scale\n";
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const Z1Chain &chain = selected[i];
        mapping << i + 1 << '\t' << chain.effective_strand_id << '\t'
                << chain.parent_molecule << '\t' << chain.parent_topology
                << '\t' << chain.status << '\t' << chain.first_node << '\t'
                << chain.second_node << '\t' << chain.first_atom << '\t'
                << chain.second_atom << '\t' << chain.effective_contour_beads
                << '\t' << chain.atoms.size() << '\t'
                << chain.reacted_site_index << '\t' << chain.next_graft_atom
                << '\t';
        for (std::size_t atom = 0; atom < chain.excluded_atoms.size(); ++atom) {
            if (atom) mapping << ',';
            mapping << chain.excluded_atoms[atom];
        }
        mapping << '\t';
        for (std::size_t atom = 0; atom < chain.atoms.size(); ++atom) {
            if (atom) mapping << ',';
            mapping << chain.atoms[atom];
        }
        mapping << '\t' << std::setprecision(15) << summary.coordinate_scale << '\n';
    }
    return summary;
}

void write_report(
    const std::filesystem::path &path, const ModelInfo &info,
    const DataFile &data, const ReducedNetwork &network,
    const std::vector<DirectionalPath> &paths, const Options &options,
    const Z1Summary &z1) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    std::map<std::string, long long> status_counts;
    std::map<std::string, long long> parent_counts;
    std::map<long long, long long> component_nodes;
    std::map<long long, long long> component_edges;
    std::map<long long, long long> degree_counts;
    std::map<std::pair<long long, long long>, long long> parallel;
    long long graph_edges = 0;
    long long self_loops = 0;
    for (const EffectiveStrand &strand : network.strands) {
        ++status_counts[strand.status];
        if (strand.status == "active") {
            ++graph_edges;
            ++component_edges[strand.graph_component];
            ++parallel[std::minmax(strand.first_node, strand.second_node)];
        } else if (strand.status == "self_loop") ++self_loops;
    }
    for (const ParentRecord &parent : network.parents) ++parent_counts[parent.state];
    for (std::size_t node = 1; node < network.nodes.size(); ++node) {
        const JunctionNode &entry = network.nodes[node];
        if (entry.component > 0) ++component_nodes[entry.component];
        ++degree_counts[entry.active_degree];
    }
    long long parallel_pairs = 0;
    for (const auto &entry : parallel) if (entry.second > 1) ++parallel_pairs;
    long long cycle_rank = self_loops;
    for (const auto &component : component_nodes)
        cycle_rank += component_edges[component.first] - component.second + 1;

    out << std::setprecision(10)
        << "PDMS generic network topology report\n"
        << "case: " << info.case_name << "\n"
        << "geometry: " << info.geometry << "\n"
        << "source: " << options.data_file << "\n"
        << "strand topology: " << info.strand_topology << "\n\n"
        << "Reduction convention\n"
        << "  linear: one original chain is one effective strand\n"
        << "  ring: cut at actually reacted sites; one reaction is a dangling loop\n"
        << "  star: each center-to-arm-end path is one effective strand\n"
        << "  grafted: consecutive reacted side-chain ends are connected by their unique molecular path\n"
        << "  crosslinkers and moderators joined by reaction bonds are collapsed into chemical junction clusters\n\n"
        << "Conversion\n"
        << "  strand sites: " << network.strand_reactive_reacted << " / "
        << network.strand_reactive_total << " ("
        << (network.strand_reactive_total == 0 ? 0.0 :
            100.0 * network.strand_reactive_reacted / network.strand_reactive_total)
        << " %)\n"
        << "  crosslinker sites: " << network.crosslinker_reactive_reacted
        << " / " << network.crosslinker_reactive_total << " ("
        << (network.crosslinker_reactive_total == 0 ? 0.0 :
            100.0 * network.crosslinker_reactive_reacted /
                network.crosslinker_reactive_total) << " %)\n\n"
        << "Reduced graph\n"
        << "  junction nodes: " << network.nodes.size() - 1 << "\n"
        << "  effective strands: " << network.strands.size() << "\n"
        << "  active graph edges: " << graph_edges << "\n"
        << "  active components: " << component_nodes.size() << "\n"
        << "  graph cycle rank: " << cycle_rank << "\n"
        << "  self-loop edges: " << self_loops << "\n"
        << "  node pairs with parallel strands: " << parallel_pairs << "\n";
    for (const auto &entry : status_counts)
        out << "  status " << entry.first << ": " << entry.second << "\n";
    out << "\nParent molecule states\n";
    for (const auto &entry : parent_counts)
        out << "  " << entry.first << ": " << entry.second << "\n";
    out << "\nActive degree distribution\n";
    for (const auto &entry : degree_counts)
        out << "  degree " << entry.first << ": " << entry.second << " nodes\n";

    out << "\nDirectional self-shortest paths to the +1 periodic image\n"
        << "  image search bound: +/-" << options.image_search_bound << " cells\n";
    const int dimensions = info.periodic_z() ? 3 : 2;
    for (int axis = 0; axis < dimensions; ++axis) {
        std::vector<double> values;
        for (const DirectionalPath &entry : paths)
            if (entry.axis == axis && std::isfinite(entry.contour_length))
                values.push_back(entry.contour_length);
        out << "  " << axis_name(axis) << ": reachable nodes " << values.size();
        if (!values.empty()) {
            const double mean = std::accumulate(values.begin(), values.end(), 0.0) /
                values.size();
            out << ", mean " << mean << " A, median " << percentile(values, 0.5)
                << " A, min " << *std::min_element(values.begin(), values.end())
                << " A, max " << *std::max_element(values.begin(), values.end()) << " A";
        }
        out << '\n';
    }
    if (!info.periodic_z())
        out << "  z omitted because the film is nonperiodic in z\n";
    if (options.skip_self_paths)
        out << "  calculation skipped by command-line request\n";

    out << "\nZ1+ export\n"
        << "  selection: " << options.z1_selection << "\n"
        << "  candidate paths selected: " << z1.selected_paths << "\n"
        << "  chains written: " << z1.chains_written << "\n"
        << "  empty trimmed strands skipped: " << z1.empty_strands_skipped << "\n"
        << "  contour bead occurrences excluded: "
        << z1.excluded_contour_occurrences << "\n"
        << "  coordinate scale: " << z1.coordinate_scale << "\n";
    if (options.z1_scaling_requested)
        out << "  scaled maximum bond target: " << options.z1_max_bond << "\n";
    else
        out << "  scaling: disabled; physical coordinates and full box lengths preserved\n";
    out << "  crosslinker geometry and reaction bonds are not part of chain contours\n"
        << "  ring reaction-site endpoints and star center beads are excluded\n"
        << "  dangling, dangling-loop, and self-loop paths are retained by the network selection\n"
        << "  grafted contours are partitioned by ordered reacted grafts without shared beads\n"
        << "  graph connectivity is retained only in the companion mapping table\n";
    if (!info.periodic_z())
        out << "  WARNING: native three-line Z1 format does not encode p p f boundaries; validate confinement/surface handling before film PPA\n";
    out << "\nBox (A): " << data.box.lx() << ' ' << data.box.ly() << ' '
        << data.box.lz() << '\n';
}

void print_help(const char *program) {
    std::cout
        << "Usage: " << program << " <case>.npt_eq <case>.info [options]\n\n"
        << "Options:\n"
        << "  --output-dir PATH          output directory (default: analysis_<case>)\n"
        << "  --z1-selection MODE        network (default), active, active-and-dangling,\n"
        << "                             or all-linearizable\n"
        << "  --z1-max-bond X            opt-in uniform scaling for PPA bond limit\n"
        << "  --image-search-bound N     lifted-cell search bound (default 2)\n"
        << "  --skip-self-paths          omit directional periodic-image searches\n"
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
        const auto value = [&]() -> std::string {
            if (++index >= argc) throw std::runtime_error("missing value for " + option);
            return argv[index];
        };
        if (option == "--output-dir") options.output_directory = value();
        else if (option == "--z1-selection") options.z1_selection = value();
        else if (option == "--z1-max-bond") {
            options.z1_max_bond = std::stod(value());
            options.z1_scaling_requested = true;
        }
        else if (option == "--image-search-bound")
            options.image_search_bound = std::stoi(value());
        else if (option == "--skip-self-paths") options.skip_self_paths = true;
        else if (option == "--help") { print_help(argv[0]); std::exit(0); }
        else throw std::runtime_error("unknown option: " + option);
    }
    if (options.z1_scaling_requested && options.z1_max_bond <= 0.0)
        throw std::runtime_error("--z1-max-bond must be positive");
    if (options.image_search_bound < 1 || options.image_search_bound > 4)
        throw std::runtime_error("--image-search-bound must be between 1 and 4");
    if (options.z1_selection != "network" &&
        options.z1_selection != "active" &&
        options.z1_selection != "active-and-dangling" &&
        options.z1_selection != "all-linearizable")
        throw std::runtime_error("invalid --z1-selection");
    return options;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        const DataFile data = parse_data_file(options.data_file, info);
        const ReducedNetwork network = reduce_network(data, info);
        const std::filesystem::path directory = analysis_directory(
            options.data_file, info, options.output_directory);
        pdms_analysis::create_directory(directory);
        const std::string name = safe_case_name(info.case_name);
        const std::vector<DirectionalPath> paths = options.skip_self_paths
            ? std::vector<DirectionalPath>{}
            : directional_self_paths(network, info.periodic_z(),
                                     options.image_search_bound);
        write_strands(directory / ("effective_strands." + name + ".tsv"), network);
        write_nodes(directory / ("network_nodes." + name + ".tsv"), network);
        write_parents(directory / ("parent_molecules." + name + ".tsv"), network);
        write_directional_paths(
            directory / ("directional_self_paths." + name + ".tsv"), paths);
        const Z1Summary z1 = write_z1(
            directory / ("config." + name + ".Z1"),
            directory / ("config." + name + ".Z1.map.tsv"),
            network, data, info, options);
        write_report(directory / ("topology_report." + name + ".txt"),
                     info, data, network, paths, options, z1);
        std::cout << "Topology analysis written to " << directory.string() << '\n';
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "topology_analyzer: " << error.what() << '\n';
        return 1;
    }
}
