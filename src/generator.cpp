#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "source_versions.hpp"
#include "composite.hpp"
#include "../vendor/Silicone_Oil/Generator/config_input.hpp"

int siliconelab_oil_main(int, char**);
int siliconelab_elastomer_main(int, char**);
int siliconelab_coating_v22_main(int, char**);
int siliconelab_coating_v35_main(int, char**);

namespace {
namespace fs = std::filesystem;
using Entry = std::pair<std::string, std::string>;
using Main = int (*)(int, char**);

struct Backend {
    Main run;
    std::string repository;
    std::string commit;
};

Backend backend(const std::string& system, const std::string& model) {
    if (system == "oil" && model.empty())
        return {siliconelab_oil_main, "https://github.com/sitengz/Silicone_Oil", oil_commit};
    if (system == "elastomer" && model.empty())
        return {siliconelab_elastomer_main, "https://github.com/sitengz/PDMS_Elastomer", network_commit};
    if (system == "coating") {
        if (model == "v22")
            return {siliconelab_coating_v22_main, "https://github.com/sitengz/Silicone_Coating", coating_commit};
        if (model == "v35")
            return {siliconelab_coating_v35_main, "https://github.com/sitengz/Silicone_Coating", coating_commit};
        throw std::runtime_error("Coating requires model = v22 or model = v35");
    }
    if (system == "oil" || system == "elastomer")
        throw std::runtime_error("model is only supported for system = coating");
    throw std::runtime_error("system must be oil, elastomer, or coating");
}

int invoke(Main run, std::vector<std::string> args) {
    std::vector<char*> argv;
    for (auto& arg : args) argv.push_back(arg.data());
    argv.push_back(nullptr);
    return run(static_cast<int>(args.size()), argv.data());
}

std::string json_quote(const std::string& value) {
    std::string result = "\"";
    constexpr char hex[] = "0123456789abcdef";
    for (unsigned char ch : value) {
        if (ch == '"' || ch == '\\') {
            result += '\\';
            result += static_cast<char>(ch);
        } else if (ch < 32) {
            result += "\\u00";
            result += hex[ch >> 4];
            result += hex[ch & 15];
        } else result += static_cast<char>(ch);
    }
    return result + '"';
}

void write_file(const fs::path& path, const std::string& text) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open output: " + path.string());
    out << text;
    out.close();
    if (!out) throw std::runtime_error("Cannot write output: " + path.string());
}

std::vector<std::string> reset_surface_momentum(const fs::path& directory,
                                               const std::string& system) {
    const std::string name = directory.filename().string();
    const fs::path input = directory / ("in." + name +
        (system == "oil" ? ".film" : ".surface"));
    // Elastomer surface inputs are built by the film workflow. Bulk coating
    // packages have no surface input; their confined preparation stays native.
    if (system == "elastomer" || !fs::is_regular_file(input)) return {};
    std::ifstream source(input);
    if (!source) throw std::runtime_error("Cannot read surface input: " + input.string());
    std::string text, line;
    bool inserted = false;
    while (std::getline(source, line)) {
        std::istringstream words(line);
        std::string command, id, group, style;
        words >> command >> id >> group >> style;
        // The first NVT follows vacuum expansion and, for oil, the brief
        // guarded lateral preparation. Do not reset between sampling phases.
        if (!inserted && command == "fix" && group == "all" && style == "nvt") {
            text += "# SiliconeLab: remove inherited COM translation before surface relaxation.\n";
            text += "velocity all zero linear\n";
            inserted = true;
        }
        text += line + '\n';
    }
    if (!source.eof() || !inserted)
        throw std::runtime_error("Cannot locate surface relaxation in " + input.string());
    write_file(input, text);
    return {input.filename().string()};
}

void write_manifest(const fs::path& directory, const fs::path& config,
                    const std::string& system, const std::string& model,
                    const Backend& selected, const std::string& data_name,
                    const std::vector<std::string>& momentum_inputs,const std::string& base) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(directory))
        if (entry.is_regular_file()) files.push_back(entry.path().filename().string());
    files.push_back("siliconelab.json");
    std::sort(files.begin(), files.end());
    std::string text = "{\n  \"schema_version\": 1,\n  \"status\": \"generated\",\n";
    text += "  \"system\": " + json_quote(system) + ",\n  \"model\": ";
    text += (model.empty() ? "null" : json_quote(model));
    text += ",\n  \"source_config\": " + json_quote(config.generic_string());
    text += ",\n  \"request_config\": \"request.conf\",\n  \"backend_config\": \"generator.conf\",";
    text += "\n  \"data_file\": " + json_quote(data_name);
    text += ",\n  \"info_file\": " + json_quote(directory.filename().string() + ".info");
    if (fs::is_regular_file(directory / "nanoparticles.json")) {
        text += ",\n  \"series\": \"composite\",\n  \"composite_file\": \"composite.json\",\n  \"nanoparticles_file\": \"nanoparticles.json\",\n  \"np_template_file\": \"np-template.csv\"";
        text += ",\n  \"base\": "+json_quote(base)+",\n  \"matrix_system\": "+json_quote(base=="elastomer"?"elastomer":"coating");
        text += ",\n  \"matrix_model\": "+(base=="elastomer"?std::string("null"):json_quote(base));
    }
    text += ",\n  \"upstream\": {\n    \"repository\": " + json_quote(selected.repository);
    text += ",\n    \"commit\": " + json_quote(selected.commit) + "\n  },\n";
    text += "  \"surface_momentum_reset\": {\"command\": \"velocity all zero linear\", ";
    text += "\"timing\": \"once before surface relaxation, after vacuum expansion and guarded preparation\", \"input_files\": [";
    for (std::size_t i = 0; i < momentum_inputs.size(); ++i)
        text += (i ? ", " : "") + json_quote(momentum_inputs[i]);
    text += "]},\n  \"generated_files\": [";
    for (std::size_t i = 0; i < files.size(); ++i)
        text += (i ? ", " : "") + json_quote(files[i]);
    text += "]\n}\n";
    write_file(directory / "siliconelab.json", text);
}

int generate(const fs::path& requested_config) {
    const fs::path config = fs::absolute(requested_config).lexically_normal();
    std::string system, model, base, output_dir;
    siliconelab_composite::Options np;
    std::set<std::string> np_keys;
    std::vector<Entry> details;
    silicone_config::read_config_file(config.string(), [&](const std::string& option, const std::string& value) {
        const auto set_once = [&](std::string& field) {
            if (!field.empty()) throw std::runtime_error("Duplicate setting: " + option);
            field = value;
        };
        if (option == "--system") set_once(system);
        else if (option == "--model") set_once(model);
        else if (option == "--base") set_once(base);
        else if (option == "--output-dir") set_once(output_dir);
        else if (option.rfind("--np-", 0) == 0) {
            if (!np_keys.insert(option).second) throw std::runtime_error("Duplicate NP setting: " + option);
            siliconelab_composite::option(np, option, value);
        }
        else if (option == "--output")
            throw std::runtime_error("Use output_dir for the common run directory instead of output");
        else if (option == "--config" || option == "--help")
            throw std::runtime_error("Reserved setting: " + option);
        else details.emplace_back(option.substr(2), value);
    });
    std::string matrix_system=system, matrix_model=model;
    if(system=="composite") {
        if(!model.empty()) throw std::runtime_error("Composite configs select base, not model");
        if(base=="v22" || base=="v35") {matrix_system="coating";matrix_model=base;}
        else if(base=="elastomer") {matrix_system="elastomer";matrix_model.clear();}
        else throw std::runtime_error("Composite requires base = v22, v35, or elastomer");
        np.supplied=true;
    } else if(!base.empty()) throw std::runtime_error("base is only supported for system = composite");
    const auto selected = backend(matrix_system, matrix_model);
    if(np.supplied) {
        if(matrix_system!="coating" && matrix_system!="elastomer")
            throw std::runtime_error("NP composites require a v22, v35, or elastomer matrix");
        siliconelab_composite::resolve(np);
        if(base.empty()) base=matrix_system=="elastomer"?"elastomer":matrix_model;
    }
    if (output_dir.empty()) throw std::runtime_error("Missing required output_dir setting");
    fs::path directory(output_dir);
    if (directory.is_relative()) directory = config.parent_path() / directory;
    directory = fs::absolute(directory).lexically_normal();
    // A trailing slash must not change the model's case name.
    if (directory.filename().empty()) directory = directory.parent_path();
    const std::string name = directory.filename().string();
    if (name.empty() || name == "." || name == "..")
        throw std::runtime_error("output_dir must name a run directory");
    // Native LAMMPS inputs interpolate case names into unquoted file tokens.
    if (name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-") != std::string::npos)
        throw std::runtime_error("Run directory name must contain only letters, digits, dots, underscores, or hyphens");
    if (directory.string().find_first_of("#\r\n") != std::string::npos)
        throw std::runtime_error("Output paths cannot contain # or line breaks (native config format)");
    const std::string data_name = "data." + name;
    // Native oil writes beside the data path; network/coating add the case name.
    const fs::path native_output = (matrix_system == "oil" ? directory : directory.parent_path()) / data_name;
    std::string native_config = "# SiliconeLab backend configuration; defaults remain upstream defaults.\n";
    for (const auto& entry : details)
        native_config += entry.first + " = \"" + entry.second + "\"\n";
    native_config += "output = \"" + native_output.generic_string() + "\"\n";

    fs::create_directories(directory.parent_path());
    if (!fs::create_directory(directory))
        throw std::runtime_error("Output directory already exists; choose a new output_dir: " + directory.string());
    try {
        fs::copy_file(config, directory / "request.conf");
        write_file(directory / "generator.conf", native_config);
        const int code = invoke(selected.run, {"siliconelab_generator", "--config", (directory / "generator.conf").string()});
        if (code != 0) throw std::runtime_error("Selected generator failed; no run package was retained");
        if (!fs::is_regular_file(directory / data_name))
            throw std::runtime_error("Selected generator did not produce the expected data file");
        if (np.supplied) siliconelab_composite::apply(directory, np, base);
        const auto momentum_inputs = reset_surface_momentum(directory, system);
        write_manifest(directory, config, system, model, selected, data_name, momentum_inputs,base);
    } catch (...) {
        // This directory was exclusively created by this call; preexisting runs
        // are rejected above and can never be removed by failure cleanup.
        std::error_code ignored;
        fs::remove_all(directory, ignored);
        throw;
    }
    std::cout << "Generated " << system << " run: " << directory.string()
              << "\nRun record: " << (directory / "siliconelab.json").string() << '\n';
    return 0;
}

void help() {
    std::cout << "Usage: siliconelab_generator --config FILE\n"
        "       siliconelab_generator --help [oil|elastomer|coating-v22|coating-v35|composite]\n"
        "       siliconelab_generator --capabilities\n\n"
        "Config format: key = value (comments start with #).\n"
        "Required common settings:\n"
        "  system = oil | elastomer | coating | composite\n"
        "  output_dir = PATH   new run directory; relative to the config file\n"
        "  model = v22 | v35   required only for coating\n"
        "  base = v22 | v35 | elastomer   required only for composite\n\n"
        "Detailed settings use the selected upstream generator's keys and defaults.\n"
        "Composites combine native matrix settings with NP settings; see --help composite.\n"
        "Existing coating/elastomer configs may also opt in with np_count.\n"
        "See --help SYSTEM for native options; use output_dir instead of output.\n"
        "Generation writes LAMMPS inputs and templates. It does not run or submit MD.\n";
}
void composite_help() {
    std::cout<<"Composite selection: system = composite; base = v22 | v35 | elastomer\n"
        "Keep the selected base's native scientific keys; output_dir must be new.\n"
        "Periodic bulk only. Representation: rigid_fcc (radius locates bead centers).\n"
        "Required: np_count (1..512)\n"
        "Geometry: np_radius=30 A; np_spacing=2^(1/6)*np_sigma; np_seed=20261014\n"
        "Beads: np_sigma=6.44584366 A; np_epsilon=1.01287845 kcal/mol; np_bead_mass=74 g/mol\n"
        "Pairs: np_interaction=wca|lj (default wca); override independently with\n"
        "  np_np_interaction, np_matrix_interaction; optional np_np_epsilon,\n"
        "  np_matrix_epsilon, np_matrix_sigma (explicit uniform matrix pair overrides).\n"
        "LJ only: np_np_cutoff, np_matrix_cutoff in A (default 2.5*sigma).\n"
        "WCA cutoff is fixed at 2^(1/6)*sigma. All NP parameters stay fixed through stages.\n"
        "Packing: np_packing_min_distance=4.5 A (intermolecular matrix beads).\n"
        "Outputs: native matrix snapshot in precursor/, composite data/input,\n"
        "  composite.json, nanoparticles.json, canonical np-template.csv, run manifest.\n"
        "Composite static analysis supports the matrix network; see docs/analyzer.md. Generation never submits jobs.\n";
}
void capabilities() {
    std::cout<<R"({"schema_version":1,"systems":["oil","elastomer","coating","composite"],"composite":{"bases":["v22","v35","elastomer"],"geometry":["bulk"],"representations":["rigid_fcc"],"pair_models":["wca","lj"],"independent_pair_channels":["np_np","np_matrix"],"elastomer_architectures":["linear","ring","star","grafted"],"coating_oil_additives":true,"np_count_limits":[1,512],"total_bead_limit":2000000,"analysis_supported":true,"analysis_scope":"bulk matrix basic/topology and composite composition/density","np_obstacle_primitive_paths_supported":false,"dynamics_analysis_supported":false},"units":{"length":"angstrom","energy":"kcal/mol","mass":"g/mol"},"side_effects":{"runs_md":false,"submits_jobs":false,"publishes":false}})"<<'\n';
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--capabilities") {capabilities();return 0;}
        if (argc == 2 && std::string(argv[1]) == "--help") { help(); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--help") {
            const std::string choice = argv[2];
            if(choice=="composite") {composite_help();return 0;}
            const auto selected = choice == "coating-v22" ? backend("coating", "v22")
                : choice == "coating-v35" ? backend("coating", "v35") : backend(choice, "");
            return invoke(selected.run, {"siliconelab_generator", "--help"});
        }
        if (argc == 3 && std::string(argv[1]) == "--config") return generate(argv[2]);
        help();
        throw std::runtime_error("Expected --config FILE or --help [SYSTEM]");
    } catch (const std::exception& error) {
        std::cerr << "SiliconeLab error: " << error.what() << '\n';
        return 1;
    }
}
