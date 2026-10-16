// Adapter for the pinned PDMS network reduction. Original vendor files are unchanged.
#pragma once
#define parse_model_info parse_upstream_model_info
#define parse_data_file parse_upstream_data_file
#include "../../vendor/PDMS_Elastomer/Analysis/network_common.hpp"
#undef parse_model_info
#undef parse_data_file
namespace pdms_analysis {
inline DataFile parse_data_file(const std::string &path, const ModelInfo &info) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open snapshot: " + path);
    std::string line;
    while (std::getline(input, line)) {
        const auto clean = trim(line);
        if (section_line(clean, "Atoms")) break;
        if (clean.find("xy xz yz") != std::string::npos ||
            clean.find("avec") != std::string::npos ||
            clean.find("bvec") != std::string::npos ||
            clean.find("cvec") != std::string::npos)
            throw std::runtime_error("triclinic snapshots are unsupported; orthorhombic bounds required");
    }
    auto data = parse_upstream_data_file(path, info);
    for (double length : {data.box.lx(), data.box.ly(), data.box.lz()})
        if (!std::isfinite(length) || length <= 0)
            throw std::runtime_error("invalid snapshot box length");
    std::set<long long> ids;
    std::set<std::pair<long long, long long>> pairs;
    for (const auto &bond : data.bonds) {
        if (bond.id < 1 || !ids.insert(bond.id).second || bond.first == bond.second ||
            !pairs.insert(std::minmax(bond.first, bond.second)).second)
            throw std::runtime_error("duplicate or invalid covalent bond");
    }
    return data;
}

inline ModelInfo parse_model_info(const std::string &path) {
    const std::string text = read_text_file(path);
    const std::string format = json_string(text, "format");
    ModelInfo info;
    info.format_version = json_integer(text, "format_version");
    if (format != "siliconelab-analysis-info" || info.format_version != 1)
        throw std::runtime_error(
            "expected normalized SiliconeLab analysis metadata version 1");
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

}
