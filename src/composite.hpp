#pragma once
// Composable bulk matrix + rigid FCC particle generator. Native sources remain unchanged.
#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include "../vendor/Silicone_Coating/Formulation/silicone_oil_component.hpp"

namespace siliconelab_composite {
namespace fs = std::filesystem;
using V = std::array<double,3>;
struct Options {
    int count = 0;
    double radius = 30;
    unsigned seed = 20261014;
    bool supplied = false;
    std::string interaction = "wca", np_interaction, matrix_interaction;
    std::string representation = "rigid_fcc";
    double sigma = silicone_oil::kDmsSigma;
    double epsilon = silicone_oil::kDmsEpsilon;
    double bead_mass = silicone_oil::kDmsMass;
    double spacing = 0, np_epsilon = 0, matrix_epsilon = 0, matrix_sigma = 0;
    double np_cutoff = 0, matrix_cutoff = 0, packing_distance = 4.5;
    double lambda = 1;
    bool native = false;
};
inline int integer(const std::string& s) {
    std::size_t end=0;
    int x=std::stoi(s,&end);
    if(end!=s.size() || x<=0) throw std::runtime_error("NP counts/seeds must be positive integers");
    return x;
}
inline double number(const std::string& s) {
    std::size_t end=0;
    double x=std::stod(s,&end);
    if(end!=s.size() || !std::isfinite(x) || x<=0) throw std::runtime_error("NP numeric settings must be positive and finite");
    return x;
}
inline void option(Options& o, const std::string& key, const std::string& value) {
    o.supplied = true;
    if(key == "--np-count") o.count=integer(value);
    else if(key == "--np-seed") o.seed=integer(value);
    else if(key == "--np-radius") o.radius=number(value);
    else if(key == "--np-representation") o.representation=value;
    else if(key == "--np-interaction") o.interaction=value;
    else if(key == "--np-np-interaction") o.np_interaction=value;
    else if(key == "--np-matrix-interaction") o.matrix_interaction=value;
    else if(key == "--np-sigma") o.sigma=number(value);
    else if(key == "--np-epsilon") o.epsilon=number(value);
    else if(key == "--np-bead-mass") o.bead_mass=number(value);
    else if(key == "--np-spacing") o.spacing=number(value);
    else if(key == "--np-np-epsilon") o.np_epsilon=number(value);
    else if(key == "--np-matrix-epsilon") o.matrix_epsilon=number(value);
    else if(key == "--np-matrix-sigma") o.matrix_sigma=number(value);
    else if(key == "--np-np-cutoff") o.np_cutoff=number(value);
    else if(key == "--np-matrix-cutoff") o.matrix_cutoff=number(value);
    else if(key == "--np-packing-min-distance") o.packing_distance=number(value);
    else if(key == "--np-lambda") o.lambda=number(value);
    else throw std::runtime_error("Unknown NP setting: "+key);
}
inline void resolve(Options& o) {
    if(o.count<1 || o.count>512) throw std::runtime_error("Composites require np_count in 1..512");
    if(o.representation!="rigid_fcc") throw std::runtime_error("np_representation currently supports rigid_fcc");
    o.native=o.interaction=="native_pdms";
    if(o.native) {
        if(!o.np_interaction.empty() || !o.matrix_interaction.empty() || o.np_epsilon || o.matrix_epsilon ||
           o.matrix_sigma || o.np_cutoff || o.matrix_cutoff || o.sigma!=silicone_oil::kDmsSigma ||
           o.epsilon!=silicone_oil::kDmsEpsilon)
            throw std::runtime_error("native_pdms inherits native stage pair parameters; use np_lambda without pair overrides");
        o.np_interaction=o.matrix_interaction="native_pdms";
    } else if(o.lambda!=1) throw std::runtime_error("np_lambda requires np_interaction = native_pdms");
    if(o.np_interaction.empty()) o.np_interaction=o.interaction;
    if(o.matrix_interaction.empty()) o.matrix_interaction=o.interaction;
    for(const auto& mode:{o.interaction,o.np_interaction,o.matrix_interaction})
        if(mode!="wca" && mode!="lj" && !(o.native && mode=="native_pdms")) throw std::runtime_error("NP interactions must be wca, lj or native_pdms");
    if((o.np_interaction=="wca" && o.np_cutoff>0) || (o.matrix_interaction=="wca" && o.matrix_cutoff>0))
        throw std::runtime_error("WCA cutoffs are derived at 2^(1/6)*sigma; explicit cutoffs require lj");
    if(o.spacing==0) o.spacing=std::pow(2.0,1.0/6.0)*o.sigma;
    if(o.radius<o.spacing) throw std::runtime_error("NP radius must include the first FCC neighbor shell");
}
struct Pair { double epsilon, sigma, cutoff; std::string interaction; };
inline Pair pair(const Options& o,int matrix_type,int np_type) {
    const bool self=matrix_type==np_type;
    const auto reference=silicone_oil::pair_parameters(1,self?1:matrix_type,300);
    Pair p{};
    p.epsilon=self ? (o.np_epsilon>0?o.np_epsilon:o.epsilon)
        : (o.matrix_epsilon>0?o.matrix_epsilon:reference.epsilon*std::sqrt(o.epsilon/silicone_oil::kDmsEpsilon));
    p.sigma=self ? o.sigma : (o.matrix_sigma>0?o.matrix_sigma:reference.sigma+0.5*(o.sigma-silicone_oil::kDmsSigma));
    p.interaction=self?o.np_interaction:o.matrix_interaction;
    if(o.native) return {reference.epsilon*o.lambda,reference.sigma,15,"native_pdms"};
    const double explicit_cut=self?o.np_cutoff:o.matrix_cutoff;
    const double minimum=std::pow(2.0,1.0/6.0)*p.sigma;
    p.cutoff=p.interaction=="wca" ? minimum : (explicit_cut>0?explicit_cut:2.5*p.sigma);
    if(!std::isfinite(p.cutoff) || !std::isfinite(p.epsilon) || !std::isfinite(p.sigma) || p.sigma<=0 || p.epsilon<=0
       || (p.interaction=="lj" && p.cutoff<=minimum))
        throw std::runtime_error("Attractive LJ requires a finite cutoff greater than 2^(1/6)*sigma");
    return p;
}
inline double maximum_cutoff(const Options& o,int np_type) {
    double cutoff=0;
    for(int i=1;i<=np_type;++i) cutoff=std::max(cutoff,pair(o,i,np_type).cutoff);
    return cutoff;
}
inline std::string read(const fs::path& p) {
    std::ifstream in(p);
    if(!in) throw std::runtime_error("Cannot read NP source: "+p.string());
    std::ostringstream out; out<<in.rdbuf();
    if(in.bad()) throw std::runtime_error("NP source read failed");
    return out.str();
}
inline void write(const fs::path& p,const std::string& s) {
    std::ofstream out(p); out<<s; out.close();
    if(!out) throw std::runtime_error("Cannot write NP output: "+p.string());
}
inline double distance2(V a,V b,double L) {
    double r2=0;
    for(int j=0;j<3;++j) {double d=a[j]-b[j]; d-=L*std::round(d/L); r2+=d*d;}
    return r2;
}
inline std::array<V,3> rotation(std::mt19937& rng) {
    std::uniform_real_distribution<double> u(0,1);
    double p=u(rng),t=2*silicone_oil::kPi*u(rng),v=2*silicone_oil::kPi*u(rng);
    double x=std::sqrt(1-p)*std::sin(t),y=std::sqrt(1-p)*std::cos(t);
    double z=std::sqrt(p)*std::sin(v),w=std::sqrt(p)*std::cos(v);
    return {{{1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)},
             {2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)},
             {2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)}}};
}
inline V rotate(V p,const std::array<V,3>& r) {
    V q{}; for(int i=0;i<3;++i) for(int j=0;j<3;++j) q[i]+=r[i][j]*p[j]; return q;
}
inline std::vector<V> fcc(double radius,double d) {
    const double step=d/std::sqrt(2.0);
    const double requested=std::ceil(radius/step);
    if(requested>64) throw std::runtime_error("NP radius exceeds local FCC construction limit");
    int n=static_cast<int>(requested);
    std::vector<V> beads;
    for(int i=-n;i<=n;++i) for(int j=-n;j<=n;++j) for(int k=-n;k<=n;++k)
        if((i+j+k)%2==0) {
            V p{step*i,step*j,step*k};
            if(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]<=radius*radius*(1+1e-12)) beads.push_back(p);
        }
    return beads;
}
struct Atom {int id,mol,type; double q; V p; std::array<int,3> image{};};
struct Grid {
    double L,cut;
    int n;
    std::unordered_map<long long,std::vector<V>> bins;
    Grid(double length,double cutoff):L(length),cut(cutoff),n(static_cast<int>(length/cutoff)) {
        if(n<3) throw std::runtime_error("Composite box too small for packing");
    }
    std::array<int,3> cell(V p) const {
        std::array<int,3> q{};
        for(int j=0;j<3;++j) q[j]=std::min(n-1,static_cast<int>((p[j]+L/2)*n/L));
        return q;
    }
    long long key(std::array<int,3> q) const {return (1LL*q[0]*n+q[1])*n+q[2];}
    bool blocked(V p) const {
        auto q=cell(p);
        for(int i=-1;i<=1;++i) for(int j=-1;j<=1;++j) for(int k=-1;k<=1;++k) {
            std::array<int,3> z{(q[0]+i+n)%n,(q[1]+j+n)%n,(q[2]+k+n)%n};
            auto it=bins.find(key(z));
            if(it!=bins.end()) for(V b:it->second) if(distance2(p,b,L)<cut*cut) return true;
        }
        return false;
    }
    void add(V p) {bins[key(cell(p))].push_back(p);}
};
inline std::string np_coefficients(const std::string& style,const Options& o,int np_type) {
    std::ostringstream out; out<<std::setprecision(17);
    out<<"# NP parameters stay fixed through all matrix temperature stages.\n";
    for(int i=1;i<=np_type;++i) {
        const auto p=pair(o,i,np_type);
        out<<"pair_coeff "<<i<<' '<<np_type<<' '<<style<<' '<<p.epsilon<<' '<<p.sigma<<' '
           <<p.cutoff<<'\n';
    }
    return out.str();
}
inline std::string input(const std::string& native,const Options& o,int np_type,const std::string& base) {
    std::istringstream in(native);
    std::ostringstream out; out<<std::setprecision(17);
    std::string line,phase="hot";
    bool velocities=false;
    bool npt=false, minimization_hold=false;
    std::set<std::string> deformation_fixes;
    const double cut=maximum_cutoff(o,np_type);
    std::string temperature_stop,damping;
    std::map<int,std::vector<std::string>> native_reference;
    while(std::getline(in,line)) {
        std::istringstream words(line);
        std::vector<std::string> w; std::string word;
        while(words>>word) w.push_back(word);
        if(w.empty()) {out<<line<<'\n'; continue;}
        if(w[0]=="pair_style") {
            native_reference.clear();
            if(o.native) {
                if(w[1]!="lj/cut" && w[1]!="lj/gromacs") throw std::runtime_error("Unsupported native PDMS pair style");
                out<<line<<'\n'; phase=w[1]=="lj/cut"?"hot":"cold";
            } else if(w[1]=="lj/cut") {
                out<<"pair_style hybrid lj/cut "<<w[2]<<" lj/cut "<<cut<<'\n'; phase="hot";
            } else if(w[1]=="lj/gromacs") {
                out<<"pair_style hybrid lj/gromacs "<<w[2]<<' '<<w[3]<<" lj/cut "<<cut<<'\n'; phase="cold";
            } else throw std::runtime_error("Unsupported native NP pair style");
        } else if(w[0]=="pair_coeff") {
            if(o.native) {
                out<<line<<'\n';
                if(w[1]=="1") native_reference[std::stoi(w[2])]=std::vector<std::string>(w.begin()+3,w.end());
            } else {
                out<<"pair_coeff "<<w[1]<<' '<<w[2]<<' '<<(phase=="hot"?"lj/cut 1":"lj/gromacs");
                for(std::size_t i=3;i<w.size();++i) out<<' '<<w[i]; out<<'\n';
            }
            if(std::stoi(w[1])==np_type-1 && w[1]==w[2]) {
                if(o.native) {
                    out<<"# NP pairs inherit this native stage, with epsilon scaled by lambda.\n";
                    for(int i=1;i<=np_type;++i) {
                        const auto& coeff=native_reference.at(i==np_type?1:i);
                        out<<"pair_coeff "<<i<<' '<<np_type<<' '<<std::stod(coeff[0])*o.lambda;
                        for(std::size_t j=1;j<coeff.size();++j) out<<' '<<coeff[j];
                        out<<'\n';
                    }
                } else {
                    out<<np_coefficients(phase=="hot"?"lj/cut 2":"lj/cut",o,np_type);
                    out<<(phase=="hot"?"pair_modify pair lj/cut 2 shift yes\n":"pair_modify pair lj/cut shift yes\n");
                }
            }
        } else if(w[0]=="comm_modify") {
            out<<"comm_modify cutoff "<<std::max(std::stod(w[2]),o.radius+2*cut)<<'\n';
        } else if(w[0]=="read_data") {
            out<<line<<"\n\ngroup nanoparticles type "<<np_type<<"\ngroup precursor subtract all nanoparticles\n"
               <<"# Rigid constituents have no internal nonbonded interactions.\n"
               <<"neigh_modify exclude molecule/intra nanoparticles\n"
               <<"compute precursor_temp precursor temp\n"
               <<"compute composite_temp all temp\n"
               <<"compute composite_press all pressure composite_temp\n";
        } else if(w[0]=="minimize") {
            out<<"# Keep the already packed rigid colloids fixed during matrix minimization.\n"
               <<"fix np_minimize_hold nanoparticles setforce 0 0 0\n"<<line<<'\n'
               <<"unfix np_minimize_hold\n";
            minimization_hold=true;
        } else if(w[0]=="velocity" && w[1]=="all" && w[2]=="create") {
            out<<"velocity precursor";
            for(std::size_t i=2;i<w.size();++i) out<<' '<<w[i]; out<<'\n';
        } else if(w[0]=="fix" && (w[1]=="integrate" || w[1]=="msd_integrate")) {
            if(!velocities) {
                out<<"# Initialize added colloids; the native precursor initialization is retained.\n"
                   <<"velocity nanoparticles create "<<w[5]<<' '<<o.seed<<" mom yes rot yes dist gaussian\n";
                velocities=true;
            }
            // LAMMPS requires a rigid fix to precede any box-changing fix.
            out<<"fix np_"<<w[1]<<" nanoparticles rigid/nvt/small molecule temp "
               <<w[5]<<' '<<w[6]<<' '<<w[7]<<" reinit no\n";
            out<<"fix "<<w[1]<<" precursor";
            for(std::size_t i=3;i<w.size();++i) out<<' '<<w[i];
            // Rigid fixes remap their centers independently. Affine rescaling
            // must affect only the flexible atoms, not NP constituent offsets.
            if(w[3]=="npt") out<<" dilate precursor";
            out<<'\n';
            if(w[3]=="npt") out<<"fix_modify "<<w[1]<<" temp precursor_temp press composite_press\n";
            npt=w[3]=="npt";
            temperature_stop=w[6]; damping=w[7];
        } else if(w[0]=="unfix" && (w[1]=="integrate" || w[1]=="msd_integrate")) {
            out<<line<<"\nunfix np_"<<w[1]<<'\n';
        } else if(w[0]=="fix" && w[1]=="xlink") {
            out<<"fix xlink precursor"; for(std::size_t i=3;i<w.size();++i) out<<' '<<w[i]; out<<'\n';
        } else if(w[0]=="fix" && w.size()>3 && w[3]=="deform") {
            deformation_fixes.insert(w[1]);
            out<<"fix "<<w[1]<<" precursor"; for(std::size_t i=3;i<w.size();++i) out<<' '<<w[i]; out<<'\n';
        } else if(w[0]=="unfix" && deformation_fixes.erase(w[1])) {
            out<<line<<'\n'
               <<"# Restore rigid coordinates/images after the last box remap.\n"
               <<"run 1 start 0 stop $(step+1:%.0f)\n";
        } else if(w[0]=="write_data" && npt) {
            // The last box remap can leave boundary-spanning constituent image
            // coordinates inconsistent with the new box. Keep the rigid state
            // and take one step at fixed volume before writing/redefining it.
            out<<"# Stop box changes and restore rigid coordinates before stage output.\n"
               <<"unfix integrate\nfix integrate precursor nvt temp "
               <<temperature_stop<<' '<<temperature_stop<<' '<<damping<<'\n'
               <<"run 1 start 0 stop $(step+1:%.0f)\n"<<line<<'\n';
            npt=false;
        } else out<<line<<'\n';
    }
    if(!velocities) throw std::runtime_error("Missing native integration stages");
    if(base=="elastomer" && !minimization_hold) throw std::runtime_error("Missing native elastomer minimization stage");
    return "# SiliconeLab composite: "+base+" matrix + rigid FCC colloids\n"+out.str();
}

inline std::string native_stage_metadata(const std::string& text,int np_type) {
    std::istringstream lines(text);std::string line,style;std::ostringstream out;
    out<<std::setprecision(17)<<'[';int stages=0,pairs=0;
    while(std::getline(lines,line)) {
        std::istringstream words(line);std::vector<std::string>w;std::string token;
        while(words>>token)w.push_back(token);
        if(w.empty())continue;
        if(w[0]=="pair_style") {
            if(stages) {if(pairs!=np_type)throw std::runtime_error("Incomplete native NP stage pairs");out<<"]},";}
            style=w[1];pairs=0;++stages;
            out<<"{\"pair_style\":\""<<style<<"\",\"global_cutoffs_A\":[";
            for(std::size_t i=2;i<w.size();++i)out<<(i==2?"":",")<<std::stod(w[i]);
            out<<"],\"pairs\":[";
        } else if(w[0]=="pair_coeff" && std::stoi(w[2])==np_type) {
            out<<(pairs++?",":"")<<"{\"type_i\":"<<std::stoi(w[1])<<",\"type_j\":"<<np_type
               <<",\"epsilon_kcal_mol\":"<<std::stod(w[3])<<",\"sigma_A\":"<<std::stod(w[4])<<",\"pair_cutoffs_A\":[";
            for(std::size_t i=5;i<w.size();++i)out<<(i==5?"":",")<<std::stod(w[i]);
            out<<"]}";
        }
    }
    if(stages!=2 || pairs!=np_type)throw std::runtime_error("Expected hot and cold native NP pair stages");
    out<<"]}]";return out.str();
}

inline void apply(const fs::path& folder,const Options& o,const std::string& base) {
    const std::string name=folder.filename().string();
    const fs::path data_path=folder/("data."+name),info_path=folder/(name+".info"),in_path=folder/("in."+name);
    const std::string data=read(data_path),native_info=read(info_path),native_input=read(in_path);
    std::string line;
    const auto marker=data.find("Atoms # full\n\n");
    const auto tail=data.find("\nBonds\n",marker);
    if(marker==std::string::npos || tail==std::string::npos) throw std::runtime_error("Unexpected native matrix data sections");
    {
        std::istringstream commands(native_input);
        bool bulk=false;
        while(std::getline(commands,line)) {
            std::istringstream row(line); std::string command,x,y,z;
            if(row>>command && command=="boundary") {row>>x>>y>>z; bulk=x=="p" && y=="p" && z=="p";}
        }
        if(!bulk) throw std::runtime_error("Composites currently require periodic bulk; positive thickness is unsupported");
    }
    const std::string head=data.substr(0,marker),topology=data.substr(tail);
    std::istringstream source_atoms(data.substr(marker+14,tail-marker-14));
    std::vector<Atom> atoms;
    std::map<int,std::vector<int>> molecules;
    int max_mol=0;
    while(std::getline(source_atoms,line)) {
        if(line.empty()) continue;
        Atom a{}; std::istringstream row(line);
        if(!(row>>a.id>>a.mol>>a.type>>a.q>>a.p[0]>>a.p[1]>>a.p[2]) || a.id!=int(atoms.size()+1) || a.type<1 || a.mol<1)
            throw std::runtime_error("Invalid native matrix atom record");
        row>>a.image[0]>>a.image[1]>>a.image[2];
        molecules[a.mol].push_back(static_cast<int>(atoms.size()));
        max_mol=std::max(max_mol,a.mol); atoms.push_back(a);
    }
    V oldbox{};
    int matrix_types=0;
    std::istringstream headers(head);
    while(std::getline(headers,line)) {
        if(line.find(" atom types")!=std::string::npos) {std::istringstream r(line);r>>matrix_types;}
        for(int j=0;j<3;++j) if(line.find(std::string(" ")+"xyz"[j]+"lo "+"xyz"[j]+"hi")!=std::string::npos) {
            std::istringstream r(line); double lo,hi; r>>lo>>hi; oldbox[j]=hi-lo;
        }
    }
    const double oldL=oldbox[0];
    if(!std::isfinite(oldL) || oldL<=0 || std::abs(oldbox[1]-oldL)>1e-6 || std::abs(oldbox[2]-oldL)>1e-6)
        throw std::runtime_error("The supported composite bases require cubic periodic bulk boxes");
    if(matrix_types!=(base=="elastomer"?3:5)) throw std::runtime_error("Unsupported native matrix type namespace");
    const int np_type=matrix_types+1;
    std::map<int,double> masses;
    const auto mass_marker=head.find("Masses\n\n");
    if(mass_marker==std::string::npos) throw std::runtime_error("Missing native mass table");
    std::istringstream mass_lines(head.substr(mass_marker+8));
    while(std::getline(mass_lines,line)) {
        std::istringstream r(line); int type; double mass;
        if(r>>type>>mass) {
            if(type<1 || type>matrix_types || !std::isfinite(mass) || mass<=0 || masses.count(type))
                throw std::runtime_error("Invalid native mass table");
            masses[type]=mass;
        }
    }
    if(int(masses.size())!=matrix_types) throw std::runtime_error("Incomplete native mass table");
    double matrix_mass=0;
    for(const auto& a:atoms) {
        if(!masses.count(a.type)) throw std::runtime_error("Atom type missing from native mass table");
        matrix_mass+=masses.at(a.type);
    }
    const double d=o.spacing;
    const auto particle=fcc(o.radius,d);
    const std::size_t base_count=atoms.size(),np_beads=particle.size()*o.count,total=base_count+np_beads;
    if(total>2000000) throw std::runtime_error("NP series exceeds the local 2 million bead limit");
    const double np_mass=np_beads*o.bead_mass;
    const double L=oldL*std::cbrt((matrix_mass+np_mass)/matrix_mass);
    double matrix_core=0;
    for(int i=1;i<=matrix_types;++i) matrix_core=std::max(matrix_core,std::pow(2.0,1.0/6.0)*pair(o,i,np_type).sigma);
    const double np_core=std::pow(2.0,1.0/6.0)*pair(o,np_type,np_type).sigma;
    const double margin=std::max(matrix_core,np_core);
    if(!std::isfinite(L) || 2*(o.radius+margin)>=L) throw std::runtime_error("NP too large for this precursor box");
    if(maximum_cutoff(o,np_type)>=L/2) throw std::runtime_error("NP pair cutoff must be smaller than half the initial box length");
    std::mt19937 rng(o.seed);
    std::uniform_real_distribution<double> center(-L/2+o.radius+margin,L/2-o.radius-margin),uniform(-L/2,L/2);
    std::vector<V> centers;
    for(int n=0;n<o.count;++n) {
        bool placed=false;
        for(int trial=0;trial<100000 && !placed;++trial) {
            V c{center(rng),center(rng),center(rng)};
            placed=true;
            for(V p:centers) if(distance2(c,p,L)<std::pow(2*o.radius+np_core,2)) {placed=false;break;}
            if(placed) centers.push_back(c);
        }
        if(!placed) throw std::runtime_error("Cannot place nonoverlapping NPs; adjust size/count/density");
    }
    // Pack whole molecules with their native internal geometry. Coordinate changes
    // cannot change native topology, reactive-site identities, or stoichiometry.
    std::vector<std::vector<int>> adjacency(base_count);
    auto bond_end=topology.size();
    for(const auto section:{"\nAngles\n","\nDihedrals\n","\nImpropers\n"}) {
        const auto pos=topology.find(section); if(pos!=std::string::npos) bond_end=std::min(bond_end,pos);
    }
    std::istringstream bond_lines(topology.substr(7,bond_end-7));
    while(std::getline(bond_lines,line)) {
        if(line.empty()) continue;
        std::istringstream r(line);int id,type,a,b;
        if(!(r>>id>>type>>a>>b) || a<1 || b<1 || a>int(base_count) || b>int(base_count) || atoms[a-1].mol!=atoms[b-1].mol)
            throw std::runtime_error("Invalid native intramolecular bond graph");
        adjacency[a-1].push_back(b-1); adjacency[b-1].push_back(a-1);
    }
    Grid grid(L,o.packing_distance);
    std::vector<V> unwrapped(base_count);
    std::vector<char> seen(base_count,0);
    for(const auto& entry:molecules) {
        const auto& ids=entry.second;
        std::vector<V> relative;
        V mean{};
        std::vector<int> queue{ids[0]}; seen[ids[0]]=1; unwrapped[ids[0]]=atoms[ids[0]].p;
        for(std::size_t q=0;q<queue.size();++q) {
            const int a=queue[q];
            for(int b:adjacency[a]) if(!seen[b]) {
                for(int j=0;j<3;++j) {
                    double delta=atoms[b].p[j]-atoms[a].p[j];
                    delta-=oldL*std::round(delta/oldL);
                    unwrapped[b][j]=unwrapped[a][j]+delta;
                }
                seen[b]=1; queue.push_back(b);
            }
        }
        if(queue.size()!=ids.size()) throw std::runtime_error("Native matrix molecule has a disconnected bond graph");
        for(int id:ids) relative.push_back(unwrapped[id]);
        for(V q:relative) for(int j=0;j<3;++j) mean[j]+=q[j]/relative.size();
        for(auto& q:relative) for(int j=0;j<3;++j) q[j]-=mean[j];
        bool placed=false;
        std::vector<V> raw,wrapped;
        for(int trial=0;trial<10000 && !placed;++trial) {
            V c{uniform(rng),uniform(rng),uniform(rng)};
            const auto rot=rotation(rng);
            raw.clear(); wrapped.clear(); placed=true;
            for(V q:relative) {
                V r=rotate(q,rot),w{};
                for(int j=0;j<3;++j) {r[j]+=c[j]; w[j]=r[j]-L*std::floor((r[j]+L/2)/L);}
                for(V np:centers) if(distance2(w,np,L)<std::pow(o.radius+matrix_core,2)) {placed=false;break;}
                if(!placed || grid.blocked(w)) {placed=false;break;}
                raw.push_back(r); wrapped.push_back(w);
            }
        }
        if(!placed) throw std::runtime_error("Cannot pack matrix molecule around NPs; adjust initial density");
        for(std::size_t i=0;i<ids.size();++i) {
            auto& a=atoms[ids[i]]; a.p=wrapped[i];
            for(int j=0;j<3;++j) a.image[j]=static_cast<int>(std::floor((raw[i][j]+L/2)/L));
            grid.add(a.p);
        }
    }
    for(int n=0;n<o.count;++n) {
        const auto rot=rotation(rng);
        for(V p:particle) {
            V r=rotate(p,rot); for(int j=0;j<3;++j) r[j]+=centers[n][j];
            atoms.push_back({static_cast<int>(atoms.size()+1),max_mol+n+1,np_type,0,r,{0,0,0}});
        }
    }
    const fs::path backup=folder/"precursor";
    fs::create_directory(backup);
    for(const auto& p:{data_path,info_path,in_path}) fs::copy_file(p,backup/p.filename());
    std::ostringstream output; output<<std::setprecision(17);
    std::string clean_head=head;
    while(!clean_head.empty() && std::isspace(static_cast<unsigned char>(clean_head.back()))) clean_head.pop_back();
    std::istringstream header_lines(clean_head);
    while(std::getline(header_lines,line)) {
        if(line.find(" atoms")!=std::string::npos) output<<total<<" atoms\n";
        else if(line.find(" atom types")!=std::string::npos) output<<np_type<<" atom types\n";
        else if(line.find("lo ")!=std::string::npos && line.find("hi")!=std::string::npos) {
            const auto axis=line.find(" xlo")!=std::string::npos?"x":line.find(" ylo")!=std::string::npos?"y":"z";
            output<<-L/2<<' '<<L/2<<' '<<axis<<"lo "<<axis<<"hi\n";
        } else output<<line<<'\n';
    }
    output<<np_type<<' '<<o.bead_mass<<" # rigid NP constituent\n\nAtoms # full\n\n";
    for(const auto& a:atoms) output<<a.id<<' '<<a.mol<<' '<<a.type<<' '<<a.q<<' '<<a.p[0]<<' '<<a.p[1]<<' '<<a.p[2]
        <<' '<<a.image[0]<<' '<<a.image[1]<<' '<<a.image[2]<<'\n';
    output<<topology;
    write(data_path,output.str());
    const auto resolved_input=input(native_input,o,np_type,base);
    write(in_path,resolved_input);
    std::ostringstream metadata; metadata<<std::setprecision(17);
    metadata<<"{\n  \"format\": \"siliconelab-composite-info\",\n  \"format_version\": 1,\n"
        <<"  \"series\": \"composite\",\n  \"base\": \""<<base<<"\",\n  \"geometry\": \"bulk\",\n"
        <<"  \"matrix_system\": \""<<(base=="elastomer"?"elastomer":"coating")<<"\",\n"
        <<"  \"representation\": \"rigid_fcc\",\n  \"np_count\": "<<o.count<<",\n"
        <<"  \"beads_per_np\": "<<particle.size()<<",\n  \"np_beads\": "<<np_beads<<",\n"
        <<"  \"precursor_beads\": "<<base_count<<",\n  \"total_beads\": "<<total<<",\n"
        <<"  \"precursor_molecules\": "<<molecules.size()<<",\n  \"total_molecules\": "<<molecules.size()+o.count<<",\n"
        <<"  \"np_molecule_id_start\": "<<max_mol+1<<",\n  \"np_atom_type\": "<<np_type<<",\n"
        <<"  \"matrix_atom_types\": "<<matrix_types<<",\n"
        <<"  \"mass_g_mol_per_bead\": "<<o.bead_mass<<",\n"
        <<"  \"matrix_mass_g_mol_equivalent\": "<<matrix_mass<<",\n"
        <<"  \"np_mass_g_mol_equivalent\": "<<np_mass<<",\n"
        <<"  \"total_mass_g_mol_equivalent\": "<<matrix_mass+np_mass<<",\n"
        <<"  \"np_weight_percent\": "<<100.0*np_mass/(matrix_mass+np_mass)<<",\n"
        <<"  \"np_radius_bead_centers_A\": "<<o.radius<<",\n  \"fcc_nearest_neighbor_A\": "<<d<<",\n"
        <<"  \"np_pair_reference_temperature_K\": 300,\n  \"sigma_A\": "<<o.sigma<<",\n"
        <<"  \"epsilon_kcal_mol\": "<<o.epsilon<<",\n  \"np_pair_cutoff_A\": "<<pair(o,np_type,np_type).cutoff<<",\n"
        <<"  \"np_np_interaction\": \""<<o.np_interaction<<"\",\n"
        <<"  \"np_matrix_interaction\": \""<<o.matrix_interaction<<"\",\n"
        <<"  \"np_lambda\": "<<o.lambda<<",\n"
        <<"  \"np_parameter_policy\": \""<<(o.native?"native_stage_matched":"fixed_300K_reference")<<"\",\n"
        <<"  \"np_native_stages\": "<<(o.native?native_stage_metadata(resolved_input,np_type):"null")<<",\n"
        <<"  \"np_pair_model\": \""<<(o.native?"native PDMS stage potential, epsilon multiplied by np_lambda":"energy-shifted constituent LJ; WCA at the minimum or attractive LJ beyond it")<<"\",\n"
        <<"  \"np_parameter_origin\": \""<<(o.native?"Exact native input coefficients; np_native_stages is authoritative":"PDMS-equivalent 300 K reference with explicitly configured overrides")<<"\",\n"
        <<"  \"np_internal_pairs\": \"excluded within each rigid molecule\",\n"
        <<"  \"packing_seed\": "<<o.seed<<",\n  \"box_length_A\": "<<L<<",\n"
        <<"  \"precursor_packing_min_distance_A\": "<<o.packing_distance<<",\n"
        <<"  \"matrix_np_center_clearance_A\": "<<o.radius+matrix_core<<",\n"
        <<"  \"np_np_center_clearance_A\": "<<2*o.radius+np_core<<",\n"
        <<"  \"placement\": \"NP-first; molecular bond-graph unwrapping, rigid rotations and periodic packing\",\n"
        <<"  \"precursor_interactions\": \"unchanged native matrix bonded terms and hot/cold pair matrices\",\n"
        <<"  \"rigid_coordinate_updates\": \"one fixed-volume step after each box-changing phase, resolved from runtime step\",\n"
        <<"  \"matrix_minimization\": \"native settings retained; NP coordinates held fixed during minimization\",\n"
        <<"  \"np_initialization\": \"added NP velocities at 800 K, seed as packing_seed\",\n"
        <<"  \"analysis_status\": \"bulk matrix basic/topology supported; composite dynamics and NP-obstacle primitive paths unsupported\",\n"
        <<"  \"precursor_info\": "<<native_info<<",\n  \"np_pairs\": [\n";
    for(int i=1;i<=np_type;++i) {
        const auto p=pair(o,i,np_type);
        metadata<<"    {\"type_i\": "<<i<<", \"type_j\": "<<np_type
                <<", \"interaction\": \""<<p.interaction<<"\", \"epsilon_kcal_mol\": "<<p.epsilon
                <<", \"sigma_A\": "<<p.sigma<<", \"cutoff_A\": "<<p.cutoff<<'}'<<(i==np_type?"\n":",\n");
    }
    metadata<<"  ],\n  \"np_centers_A\": [\n";
    for(std::size_t i=0;i<centers.size();++i) metadata<<"    ["<<centers[i][0]<<", "<<centers[i][1]<<", "<<centers[i][2]<<']'<<(i+1==centers.size()?"\n":",\n");
    metadata<<"  ]\n}\n";
    write(info_path,metadata.str());
    write(folder/"nanoparticles.json",metadata.str());
    write(folder/"composite.json",metadata.str());
    std::ostringstream coordinates;coordinates<<std::setprecision(17)<<"bead,x_A,y_A,z_A\n";
    for(std::size_t i=0;i<particle.size();++i) coordinates<<i+1<<','<<particle[i][0]<<','<<particle[i][1]<<','<<particle[i][2]<<'\n';
    write(folder/"np-template.csv",coordinates.str());
}
} // namespace siliconelab_composite
