#include "common.hpp"
#define main upstream_basic_main
#include "../../vendor/PDMS_Elastomer/Analysis/basic_network_analyzer.cpp"
#undef main

// Geometry and graph calculations reuse the pinned, independently tested reduction.
// The new reader validates each component, including variable-length oil molecules.
int main(int argc, char **argv) {
    try {
        if (argc != 4) throw std::runtime_error("expected DATA NORMALIZED_INFO OUTPUT_DIR");
        const ModelInfo info = parse_model_info(argv[2]);
        const std::string metadata = read_text_file(argv[2]);
        const DataFile data = parse_data_file(argv[1], info);
        const bool network_enabled = json_integer(metadata, "network_enabled") != 0;
        const auto expected_counts = json_integer_array(metadata, "expected_component_beads");
        const auto lower_sizes = json_integer_array(metadata, "minimum_molecule_beads");
        const auto upper_sizes = json_integer_array(metadata, "maximum_molecule_beads");
        if (expected_counts.size()!=4 || lower_sizes.size()!=4 || upper_sizes.size()!=4)
            throw std::runtime_error("invalid component validation metadata");
        if (data.declared_atoms != info.total_beads)
            throw std::runtime_error("snapshot atom count differs from metadata");
        const auto counts = component_beads(data, info);
        for (int c=0; c<4; ++c) {
            if (counts[c] != expected_counts[c])
                throw std::runtime_error("component bead count differs from metadata");
            const auto &component = info.components[c];
            for (long long m=component.molecule_start; m>0 && m<=component.molecule_end; ++m) {
                const auto n = static_cast<long long>(data.molecule_atoms[m].size());
                if (n<lower_sizes[c] || n>upper_sizes[c])
                    throw std::runtime_error("molecule bead count differs from metadata");
            }
        }
        // Reaction bonds are exclusively intermolecular crosslinker-to-functional-site bonds.
        long long reactions=0;
        std::set<long long> reacted_sites;
        for (const auto &bond:data.bonds) {
            const auto &a=data.atoms[bond.first], &b=data.atoms[bond.second];
            if (bond.type!=info.crosslink_bond_type) {
                if (a.molecule!=b.molecule)
                    throw std::runtime_error("unexpected intermolecular nonreaction bond");
                continue;
            }
            ++reactions;
            const int ca=component_for_molecule(a.molecule,info);
            const int cb=component_for_molecule(b.molecule,info);
            const bool valid=(ca==kCrosslinker && (cb==kStrand || cb==kModerator)) ||
                             (cb==kCrosslinker && (ca==kStrand || ca==kModerator));
            if (!network_enabled || !valid || a.molecule==b.molecule || a.type!=1 || b.type!=1)
                throw std::runtime_error("invalid reaction bond or reacted atom type");
            if (!reacted_sites.insert(a.id).second || !reacted_sites.insert(b.id).second)
                throw std::runtime_error("functional site reacted more than once");
            for (const Atom *atom:{&a,&b}) {
                const int c=component_for_molecule(atom->molecule,info);
                const auto rank=local_rank(atom->id,data.molecule_atoms[atom->molecule]);
                const auto &sites=c==kCrosslinker?info.crosslinker_reactive_bead_sites:info.reactive_bead_sites;
                if (c!=kModerator && std::find(sites.begin(),sites.end(),rank)==sites.end())
                    throw std::runtime_error("reaction at undeclared functional site");
                if (c==kModerator && rank==1)
                    throw std::runtime_error("moderator center cannot react");
            }
        }
        if (data.declared_bonds != json_integer(metadata,"initial_bonds")+reactions)
            throw std::runtime_error("final bond count differs from initial bonds plus reactions");
        long long unreacted2=0,unreacted3=0;
        for (const auto &a:data.atoms) if (a.seen) {
            if (a.type==2) ++unreacted2;
            if (a.type==3) ++unreacted3;
        }
        if (network_enabled && (unreacted2+reactions!=json_integer(metadata,"initial_type2_sites") ||
            unreacted3+reactions!=json_integer(metadata,"initial_type3_sites")))
            throw std::runtime_error("reactive-site conservation check failed");
        std::array<long long,6> types{{0,0,0,0,0,0}};
        for(const auto &a:data.atoms) if(a.seen) {
            if(a.type<1 || a.type>5) throw std::runtime_error("unsupported bead type");
            ++types[a.type];
        }
        if(types[4]!=json_integer(metadata,"expected_mps_backbone") ||
           types[5]!=json_integer(metadata,"expected_mps_pendant"))
            throw std::runtime_error("MPS backbone/pendant counts differ from metadata");
        const auto expected_masses=json_object_for_key(metadata,"expected_type_masses");
        for(std::size_t type=1;type<data.masses.size();++type) {
            if(data.masses[type]==0) continue;
            if(std::abs(data.masses[type]-json_number(expected_masses,std::to_string(type)))>1e-5)
                throw std::runtime_error("bead mass differs from info-file force field");
        }

        const std::filesystem::path dir(argv[3]);
        pdms_analysis::create_directory(dir);
        const auto name=safe_case_name(info.case_name);
        std::ofstream stats(dir/"snapshot_statistics.tsv");
        stats<<"quantity\tvalue\tunit\n"<<std::setprecision(15);
        const double volume=data.box.lx()*data.box.ly()*data.box.lz();
        stats<<"atoms\t"<<data.declared_atoms<<"\tcount\n"
             <<"molecules\t"<<info.total_molecules<<"\tcount\n"
             <<"volume\t"<<volume<<"\tA^3\n"
             <<"density\t"<<data.total_mass_g_per_mol/(0.602214076*volume)<<"\tg/cm^3\n"
             <<"velocity_temperature\t"<<velocity_temperature(data)<<"\tK\n"
             <<"reaction_bonds\t"<<reactions<<"\tcount\n"
             <<"box_Lx\t"<<data.box.lx()<<"\tA\n"
             <<"box_Ly\t"<<data.box.ly()<<"\tA\n"
             <<"box_Lz\t"<<data.box.lz()<<"\tA\n";
        const auto masses=component_masses(data,info);
        const std::array<std::string,4> labels{{"network","crosslinker","moderator","oil"}};
        for(int c=0;c<4;++c) {
            stats<<labels[c]<<"_beads\t"<<counts[c]<<"\tcount\n"
                 <<labels[c]<<"_mass_fraction\t"<<masses[c]/data.total_mass_g_per_mol<<"\tfraction\n";
        }
        stats.close();
        if (network_enabled) {
            DisjointSet sets(static_cast<std::size_t>(info.total_molecules+1));
            std::vector<long long> degrees(static_cast<std::size_t>(info.total_molecules+1),0);
            std::array<long long,4> reacted{{0,0,0,0}};
            for(const auto &bond:data.bonds) if(bond.type==info.crosslink_bond_type) {
                const auto &a=data.atoms[bond.first],&b=data.atoms[bond.second];
                sets.unite(a.molecule,b.molecule);++degrees[a.molecule];++degrees[b.molecule];
                ++reacted[component_for_molecule(a.molecule,info)];
                ++reacted[component_for_molecule(b.molecule,info)];
            }
            std::map<long long,std::array<long long,4>> component_molecules;
            std::map<long long,double> component_mass;
            std::map<long long,long long> component_atoms;
            std::ofstream molecule_degrees(dir/"reactive_molecule_degrees.tsv");
            molecule_degrees<<"molecule\tcomponent\treaction_degree\tchemical_component\n";
            double reactive_mass=0;long long reactive_beads=0;
            for(long long m=1;m<=info.total_molecules;++m) {
                const int c=component_for_molecule(m,info);
                if(c==kFiller) continue;
                const auto root=sets.find(m);
                ++component_molecules[root][c];
                molecule_degrees<<m<<'\t'<<labels[c]<<'\t'<<degrees[m]<<'\t'<<root<<'\n';
                for(auto a:data.molecule_atoms[m]) {
                    component_mass[root]+=data.atoms[a].mass;
                    ++component_atoms[root];reactive_mass+=data.atoms[a].mass;++reactive_beads;
                }
            }
            long long largest=0,linked=0,isolated=0;
            for(const auto &entry:component_molecules) {
                const auto n=std::accumulate(entry.second.begin(),entry.second.end(),0LL);
                if(n>1) ++linked;else ++isolated;
                if(largest==0 || component_atoms[entry.first]>component_atoms[largest]) largest=entry.first;
            }
            std::ofstream chemistry(dir/"chemical_network_statistics.tsv");
            chemistry<<"quantity\tvalue\tunit\n"<<std::setprecision(15);
            const auto initial2=json_integer(metadata,"initial_type2_sites"),initial3=json_integer(metadata,"initial_type3_sites");
            const auto strand_total=info.components[kStrand].molecules*info.strand_functionality;
            const auto moderator_total=initial2-strand_total;
            const auto fraction=[](double a,double b){return b>0?a/b:std::numeric_limits<double>::quiet_NaN();};
            chemistry<<"overall_conversion\t"<<fraction(reactions,std::min(initial2,initial3))<<"\tfraction\n"
                <<"strand_site_conversion\t"<<fraction(reacted[kStrand],strand_total)<<"\tfraction\n"
                <<"crosslinker_site_conversion\t"<<fraction(reacted[kCrosslinker],initial3)<<"\tfraction\n"
                <<"moderator_site_conversion\t"<<fraction(reacted[kModerator],moderator_total)<<"\tfraction\n"
                <<"linked_chemical_components\t"<<linked<<"\tcount\n"
                <<"isolated_reactive_molecules\t"<<isolated<<"\tcount\n"
                <<"largest_chemical_component_bead_fraction\t"<<fraction(component_atoms[largest],reactive_beads)<<"\tfraction\n"
                <<"largest_chemical_component_mass_fraction\t"<<fraction(component_mass[largest],reactive_mass)<<"\tfraction\n";
            const auto network=reduce_network(data,info);
            const auto properties=strand_properties(network,data,info);
            const auto graph=graph_summary(network);
            write_strand_properties(dir/("strand_properties."+name+".tsv"),properties);
            write_statistics(dir/("strand_statistics."+name+".tsv"),properties,graph,info);
            write_junctions(dir/("junction_properties."+name+".tsv"),network);
            write_histograms(dir/("strand_histograms."+name+".tsv"),properties,graph,40);
            write_network_statistics(dir/("network_statistics."+name+".tsv"),
                data,info,network,graph,properties,info.final_temperature_k);
            Options options;
            write_report(dir/("basic_network_report."+name+".txt"),options,data,info,
                network,graph,properties,info.final_temperature_k);
        }
        // Reconstruct each oil molecule using its covalent graph, not image flags
        // or distance from the first bead. Rg includes MPS pendants; Lc/Ree and Z1
        // use backbone beads only. Molecule shape is bead-number weighted.
        std::vector<std::vector<long long>> adjacency(data.atoms.size());
        for (const auto &bond:data.bonds) {
            if (data.atoms[bond.first].molecule!=data.atoms[bond.second].molecule) continue;
            adjacency[bond.first].push_back(bond.second);
            adjacency[bond.second].push_back(bond.first);
        }
        std::ofstream shape(dir/"oil_chain_properties.tsv"),coords(dir/"oil_backbones.tsv");
        shape<<"molecule\tbeads\tbackbone_beads\tmass_g_mol\tLc_A\tRee_x_A\tRee_y_A\tRee_z_A\tRee_A\tRg_A\tgxx_A2\tgyy_A2\tgzz_A2\tgxy_A2\tgxz_A2\tgyz_A2\tlambda1_A2\tlambda2_A2\tlambda3_A2\tanisotropy\tstraightness\ttortuosity\tP2_x\tP2_y\tP2_z\n";
        coords<<"molecule\tatom\tx_A\ty_A\tz_A\n";
        shape<<std::setprecision(15);coords<<std::setprecision(15);
        const auto &oil=info.components[kFiller];
        for(long long m=oil.molecule_start;m>0 && m<=oil.molecule_end;++m) {
            const auto &atoms=data.molecule_atoms[m];
            std::map<long long,Vec3> positions;
            std::queue<long long> q;
            positions[atoms.front()]=data.atoms[atoms.front()].position;q.push(atoms.front());
            while(!q.empty()) {
                const auto a=q.front();q.pop();
                for(auto b:adjacency[a]) if(!positions.count(b)) {
                    positions[b]=positions[a]+minimum_image_vector(
                        data.atoms[b].position-data.atoms[a].position,data.box,info.periodic_z());
                    q.push(b);
                }
            }
            if(positions.size()!=atoms.size()) throw std::runtime_error("disconnected oil molecule");
            std::vector<long long> backbone;
            std::map<long long,std::vector<long long>> bg;
            std::vector<Vec3> all;
            double mass=0;
            for(auto a:atoms) {
                all.push_back(positions[a]);mass+=data.atoms[a].mass;
                if(data.atoms[a].type!=5) {
                    backbone.push_back(a);
                    for(auto b:adjacency[a]) if(data.atoms[b].type!=5) bg[a].push_back(b);
                } else if(adjacency[a].size()!=1 || data.atoms[adjacency[a][0]].type!=4)
                    throw std::runtime_error("invalid MPS pendant connectivity");
            }
            std::vector<long long> ordered;
            if(backbone.size()==1) ordered=backbone;
            else {
                std::vector<long long> ends;
                for(auto a:backbone) {
                    if(bg[a].size()==1) ends.push_back(a);
                    if(bg[a].size()>2) throw std::runtime_error("branched oil backbone");
                }
                if(ends.size()!=2) throw std::runtime_error("oil backbone must be an open linear contour");
                long long prev=0,cur=ends.front();
                while(cur) {
                    ordered.push_back(cur);long long next=0;
                    for(auto b:bg[cur]) if(b!=prev) next=b;
                    prev=cur;cur=next;
                    if(ordered.size()>backbone.size()) throw std::runtime_error("cyclic oil contour");
                }
                if(ordered.size()!=backbone.size()) throw std::runtime_error("disconnected oil backbone");
            }
            double lc=0;
            for(std::size_t i=0;i<ordered.size();++i) {
                const auto p=positions[ordered[i]];
                coords<<m<<'\t'<<ordered[i]<<'\t'<<p.x<<'\t'<<p.y<<'\t'<<p.z<<'\n';
                if(i) lc+=norm(p-positions[ordered[i-1]]);
            }
            const auto re=positions[ordered.back()]-positions[ordered.front()];
            const double ree=norm(re),nan=std::numeric_limits<double>::quiet_NaN();
            const auto s=shape_data(all);
            shape<<m<<'\t'<<atoms.size()<<'\t'<<ordered.size()<<'\t'<<mass<<'\t'<<lc
                 <<'\t'<<re.x<<'\t'<<re.y<<'\t'<<re.z<<'\t'<<ree<<'\t'<<s.rg
                 <<'\t'<<s.gxx<<'\t'<<s.gyy<<'\t'<<s.gzz<<'\t'<<s.gxy<<'\t'<<s.gxz<<'\t'<<s.gyz
                 <<'\t'<<s.eigenvalues[0]<<'\t'<<s.eigenvalues[1]<<'\t'<<s.eigenvalues[2]
                 <<'\t'<<s.relative_shape_anisotropy<<'\t'<<(lc>0?ree/lc:nan)<<'\t'<<(ree>0?lc/ree:nan)
                 <<'\t'<<(ree>0?0.5*(3*re.x*re.x/(ree*ree)-1):nan)
                 <<'\t'<<(ree>0?0.5*(3*re.y*re.y/(ree*ree)-1):nan)
                 <<'\t'<<(ree>0?0.5*(3*re.z*re.z/(ree*ree)-1):nan)<<'\n';
        }
        if(!shape || !coords) throw std::runtime_error("failed writing molecular properties");
        std::cout<<"Validated snapshot; static structural analysis written to "<<dir.string()<<'\n';
        return 0;
    } catch(const std::exception &e) {
        std::cerr<<"siliconelab_static: "<<e.what()<<'\n';return 1;
    }
}
