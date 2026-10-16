// Export the pinned native reduction's exact contour ownership for spatial analysis.
#include "common.hpp"
#include <iostream>

int main(int argc, char **argv) {
    try {
        if (argc != 4) throw std::runtime_error("usage: siliconelab_interface_network MATRIX_DATA NORMALIZED_INFO NEW_OUTPUT_DIR");
        using namespace pdms_analysis;
        const auto info = parse_model_info(argv[2]);
        const auto data = parse_data_file(argv[1], info);
        const auto network = reduce_network(data, info);
        const std::filesystem::path out(argv[3]);
        if (!std::filesystem::create_directory(out)) throw std::runtime_error("choose a new output directory");
        std::ofstream strands(out / "contours.tsv"), nodes(out / "nodes.tsv");
        strands << "strand_id\tparent_molecule\tstatus\tfirst_node\tsecond_node\tatom_ids\n";
        for (const auto &s : network.strands) {
            strands << s.id << '\t' << s.parent_molecule << '\t' << s.status << '\t' << s.first_node << '\t' << s.second_node << '\t';
            for (std::size_t i = 0; i < s.atoms.size(); ++i) strands << (i ? "," : "") << s.atoms[i];
            strands << '\n';
        }
        nodes << std::setprecision(17) << "node_id\tkind\tactive_degree\tdangling_degree\tx_A\ty_A\tz_A\tmember_molecules\tattachment_atoms\n";
        for (std::size_t i = 1; i < network.nodes.size(); ++i) {
            const auto &n = network.nodes[i];
            nodes << n.id << '\t' << n.kind << '\t' << n.active_degree << '\t' << n.dangling_degree << '\t' << n.position.x << '\t' << n.position.y << '\t' << n.position.z << '\t';
            for (std::size_t j = 0; j < n.member_molecules.size(); ++j) nodes << (j ? "," : "") << n.member_molecules[j];
            nodes << '\t';
            for (std::size_t j = 0; j < n.attachment_atoms.size(); ++j) nodes << (j ? "," : "") << n.attachment_atoms[j];
            nodes << '\n';
        }
        if (!strands || !nodes) throw std::runtime_error("failed to write contour export");
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "interface_network: " << e.what() << '\n';
        return 1;
    }
}
