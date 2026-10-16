// General component profiles and sampled-origin dynamics; normalized metadata
// and snapshot validation use the same pinned parser as part 1.
#include "common.hpp"
#include <deque>
#include <iostream>
#include <tuple>
using namespace pdms_analysis;
namespace {
const char *labels[] = {"network", "crosslinker", "moderator", "oil"};
struct Point {
  double x = 0, y = 0, z = 0;
};
Point sub(Point a, Point b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
struct Frame {
  long long step = 0;
  Box box;
  std::vector<Point> positions;
  std::vector<long long> molecules;
  std::vector<int> types;
};
void need(bool yes, const std::string &s) {
  if (!yes)
    throw std::runtime_error(s);
}
Point centroid(const Frame &f, const DataFile &data) {
  Point s;
  double m = 0;
  for (size_t i = 1; i < f.positions.size(); ++i) {
    double w = data.atoms[i].mass;
    m += w;
    s.x += w * f.positions[i].x;
    s.y += w * f.positions[i].y;
    s.z += w * f.positions[i].z;
  }
  return {s.x / m, s.y / m, s.z / m};
}
bool frame(std::istream &in, Frame &f, const DataFile &data,
           const ModelInfo &info) {
  std::string line;
  if (!std::getline(in, line))
    return false;
  need(line == "ITEM: TIMESTEP", "trajectory frame header");
  std::getline(in, line);
  f.step = std::stoll(line);
  std::getline(in, line);
  need(line == "ITEM: NUMBER OF ATOMS", "trajectory count header");
  std::getline(in, line);
  need(std::stoll(line) == info.total_beads, "trajectory atom count mismatch");
  std::getline(in, line);
  need(line == std::string("ITEM: BOX BOUNDS pp pp ") +
                   (info.periodic_z() ? "pp" : "ff"),
       "trajectory boundary mismatch or unsupported triclinic box");
  double *lows[] = {&f.box.xlo, &f.box.ylo, &f.box.zlo};
  double *highs[] = {&f.box.xhi, &f.box.yhi, &f.box.zhi};
  for (int a = 0; a < 3; ++a) {
    std::getline(in, line);
    std::istringstream row(line);
    need(bool(row >> *lows[a] >> *highs[a]), "invalid bounds");
    need(*highs[a] > *lows[a], "invalid box lengths");
  }
  std::getline(in, line);
  std::istringstream head(line);
  std::vector<std::string> columns;
  std::string v;
  head >> v >> v;
  while (head >> v)
    columns.push_back(v);
  auto col = [&](std::string s) {
    auto p = std::find(columns.begin(), columns.end(), s);
    need(p != columns.end(), "missing trajectory column " + s);
    return int(p - columns.begin());
  };
  int id = col("id"), mol = col("mol"), type = col("type");
  bool unwrap =
      std::find(columns.begin(), columns.end(), "xu") != columns.end();
  int x = col(unwrap ? "xu" : "x"), y = col(unwrap ? "yu" : "y"),
      z = col(unwrap ? "zu" : "z"), ix = 0, iy = 0, iz = 0;
  if (!unwrap) {
    ix = col("ix");
    iy = col("iy");
    iz = col("iz");
  }
  f.positions.assign(info.total_beads + 1, {});
  f.molecules.assign(info.total_beads + 1, 0);
  f.types.assign(info.total_beads + 1, 0);
  std::vector<bool> seen(info.total_beads + 1, false);
  for (long long n = 0; n < info.total_beads; ++n) {
    need(bool(std::getline(in, line)), "truncated trajectory");
    std::istringstream row(line);
    std::vector<double> vals;
    double q;
    while (row >> q)
      vals.push_back(q);
    need(vals.size() == columns.size(), "trajectory column count");
    for (double a : vals)
      need(std::isfinite(a), "nonfinite trajectory coordinate");
    long long i = std::llround(vals[id]);
    need(i >= 1 && i <= info.total_beads && !seen[i] && vals[id] == i,
         "duplicate/invalid atom ID");
    seen[i] = true;
    f.molecules[i] = std::llround(vals[mol]);
    f.types[i] = std::llround(vals[type]);
    need(vals[mol] == f.molecules[i] && vals[type] == f.types[i],
         "noninteger trajectory identity");
    need(f.molecules[i] == data.atoms[i].molecule &&
             f.types[i] == data.atoms[i].type,
         "trajectory identity changed");
    f.positions[i] = {vals[x], vals[y], vals[z]};
    if (!unwrap) {
      need(vals[ix] == std::llround(vals[ix]) &&
               vals[iy] == std::llround(vals[iy]) &&
               vals[iz] == std::llround(vals[iz]),
           "noninteger image flag");
      f.positions[i].x += vals[ix] * f.box.lx();
      f.positions[i].y += vals[iy] * f.box.ly();
      f.positions[i].z += vals[iz] * f.box.lz();
    }
  }
  return true;
}
struct Acc {
  long long n = 0;
  double x = 0, y = 0, z = 0;
  void add(Point p) {
    ++n;
    x += p.x * p.x;
    y += p.y * p.y;
    z += p.z * p.z;
  }
};
struct Origin {
  Frame f;
  Point center;
  std::vector<Point> com;
  std::vector<int> layers;
  int first_layer = 0, layer_count = 0;
};
std::vector<Point> molecules(const Frame &f, const DataFile &d,
                             const ModelInfo &info) {
  std::vector<Point> p(info.total_molecules + 1);
  std::vector<double> masses(p.size());
  for (size_t i = 1; i < f.positions.size(); ++i) {
    auto m = f.molecules[i];
    auto a = f.positions[i];
    double w = d.atoms[i].mass;
    p[m].x += a.x * w;
    p[m].y += a.y * w;
    p[m].z += a.z * w;
    masses[m] += w;
  }
  for (size_t m = 1; m < p.size(); ++m) {
    need(masses[m] > 0, "empty molecule");
    p[m].x /= masses[m];
    p[m].y /= masses[m];
    p[m].z /= masses[m];
  }
  return p;
}
void dynamics(const std::string &path, const ModelInfo &info, const DataFile &d,
              const std::filesystem::path &out, int stride, int maxorig) {
  std::ifstream in(path);
  need(bool(in), "cannot read trajectory");
  Frame f;
  std::vector<Origin> origins;
  std::map<std::tuple<long long, int, int, int>, Acc> sums;
  long long previous = -1;
  int frames = 0;
  Box reference;
  while (frame(in, f, d, info)) {
    need(f.step > previous, "trajectory timesteps must increase; phases cannot "
                            "be concatenated across reset");
    previous = f.step;
    if (!frames)
      reference = f.box;
    need(std::abs(reference.lx() - f.box.lx()) < 1e-6 &&
             std::abs(reference.ly() - f.box.ly()) < 1e-6 &&
             std::abs(reference.lz() - f.box.lz()) < 1e-6,
         "MSD requires constant production box");
    auto center = centroid(f, d);
    auto com = molecules(f, d, info);
    if (frames % stride == 0 && int(origins.size()) < maxorig) {
      Origin origin{f, center, com, {}};
      origin.layers.resize(f.positions.size());
      int minimum = 1000000, maximum = -1000000;
      for (size_t i = 1; i < f.positions.size(); ++i) {
        int layer = int(std::floor((f.positions[i].z - center.z) / 5));
        origin.layers[i] = layer;
        minimum = std::min(minimum, layer);
        maximum = std::max(maximum, layer);
      }
      origin.first_layer = minimum;
      origin.layer_count = maximum - minimum + 1;
      need(origin.layer_count < 100000, "unreasonable origin layer extent");
      origins.push_back(std::move(origin));
    }
    for (const auto &o : origins) {
      long long lag = f.step - o.f.step;
      Point drift = sub(center, o.center);
      std::vector<std::array<Acc, 4>> layers(o.layer_count);
      std::array<Acc, 4> beads{}, molecule{};
      for (size_t i = 1; i < f.positions.size(); ++i) {
        int c = component_for_molecule(f.molecules[i], info);
        Point delta = sub(sub(f.positions[i], o.f.positions[i]), drift);
        layers[o.layers[i] - o.first_layer][c].add(delta);
        beads[c].add(delta);
      }
      for (size_t m = 1; m < com.size(); ++m) {
        int c = component_for_molecule(m, info);
        molecule[c].add(sub(sub(com[m], o.com[m]), drift));
      }
      auto merge = [&](int c, int p, int layer, const Acc &a) {
        if (!a.n)
          return;
        auto &target = sums[{lag, c, p, layer}];
        target.n += a.n;
        target.x += a.x;
        target.y += a.y;
        target.z += a.z;
      };
      for (int c = 0; c < 4; ++c) {
        merge(c, 0, 1000000, beads[c]);
        merge(c, 1, 1000000, molecule[c]);
        for (int l = 0; l < o.layer_count; ++l)
          merge(c, 0, l + o.first_layer, layers[l][c]);
      }
    }
    ++frames;
  }
  need(frames >= 2, "insufficient trajectory frames");
  std::ofstream csv(out / "dynamics.tsv");
  csv << "lag_steps\tlag_ps\tcomponent\tparticle\tlayer_center_relative_"
         "A\tobservations\tMSD_x_A2\tMSD_y_A2\tMSD_z_A2\tMSD_xy_A2\tMSD_3d_A2\n"
      << std::setprecision(14);
  for (const auto &entry : sums) {
    auto [lag, c, p, layer] = entry.first;
    const auto &a = entry.second;
    csv << lag << '\t' << lag * info.timestep_fs * .001 << '\t' << labels[c]
        << '\t' << (p ? "molecular_COM" : "beads") << '\t';
    if (layer == 1000000)
      csv << "global";
    else
      csv << (layer + .5) * 5;
    csv << '\t' << a.n << '\t' << a.x / a.n << '\t' << a.y / a.n << '\t'
        << a.z / a.n << '\t' << (a.x + a.y) / a.n << '\t'
        << (a.x + a.y + a.z) / a.n << '\n';
  }
  std::ofstream report(out / "sampling.json");
  report << "{\"frames\":" << frames
         << ",\"first_step\":" << origins.front().f.step
         << ",\"last_step\":" << previous << ",\"origins\":" << origins.size()
         << ",\"origin_stride_frames\":" << stride
         << ",\"timestep_fs\":" << info.timestep_fs
         << ",\"drift_removal\":\"whole-system mass-weighted "
            "COM\",\"averaging\":\"sampled origins pooled by actual timestep "
            "lag; layers assigned at each origin relative to material COM\"}\n";
}
void spatial(const ModelInfo &info, const DataFile &d,
             const std::filesystem::path &out) {
  int bins = std::max(1, int(std::ceil(d.box.lz() / 2.0)));
  double width = d.box.lz() / bins;
  std::vector<std::array<double, 4>> mass(bins);
  std::vector<std::array<long long, 4>> count(bins);
  for (const auto &a : d.atoms)
    if (a.seen) {
      double z = a.position.z;
      if (info.periodic_z())
        z -= std::floor((z - d.box.zlo) / d.box.lz()) * d.box.lz();
      int b = std::clamp(int((z - d.box.zlo) / width), 0, bins - 1);
      int c = component_for_molecule(a.molecule, info);
      mass[b][c] += a.mass;
      ++count[b][c];
    }
  std::ofstream f(out / "component_profiles.tsv");
  f << "z_A\tbin_width_A\tcomponent\tbeads\tmass_g_per_mol\tdensity_g_cm3\n"
    << std::setprecision(14);
  for (int b = 0; b < bins; ++b)
    for (int c = 0; c < 4; ++c)
      f << d.box.zlo + (b + .5) * width << '\t' << width << '\t' << labels[c]
        << '\t' << count[b][c] << '\t' << mass[b][c] << '\t'
        << mass[b][c] / (.602214076 * d.box.lx() * d.box.ly() * width) << '\n';
  // Export precise wrapped snapshot coordinates for chemistry/component fields.
  std::ofstream atoms(out / "snapshot_atoms.tsv");
  atoms << "id\tmolecule\tcomponent\ttype\tmass\tx\ty\tz\n"
        << std::setprecision(14);
  for (const auto &a : d.atoms)
    if (a.seen)
      atoms << a.id << '\t' << a.molecule << '\t'
            << labels[component_for_molecule(a.molecule, info)] << '\t'
            << a.type << '\t' << a.mass << '\t' << a.position.x << '\t'
            << a.position.y << '\t' << a.position.z << '\n';
  std::ofstream b(out / "box.json");
  b << std::setprecision(15) << "{\"lo\":[" << d.box.xlo << ',' << d.box.ylo
    << ',' << d.box.zlo << "],\"hi\":[" << d.box.xhi << ',' << d.box.yhi << ','
    << d.box.zhi << "],\"geometry\":\"" << info.geometry << "\"}\n";
}
} // namespace
int main(int argc, char **argv) {
  try {
    need(argc == 5 || argc == 8,
         "usage: BACKEND spatial DATA INFO OUT | BACKEND dynamics DATA INFO "
         "OUT TRAJECTORY ORIGIN_STRIDE MAX_ORIGINS");
    std::string mode = argv[1];
    auto info = parse_model_info(argv[3]);
    auto data = parse_data_file(argv[2], info);
    need(data.declared_atoms == info.total_beads,
         "snapshot/metadata count mismatch");
    std::filesystem::path out = argv[4];
    std::filesystem::create_directories(out);
    if (mode == "spatial" && argc == 5)
      spatial(info, data, out);
    else if (mode == "dynamics" && argc == 8) {
      int stride = std::stoi(argv[6]), maxorig = std::stoi(argv[7]);
      need(stride > 0 && maxorig > 0, "positive origin sampling required");
      dynamics(argv[5], info, data, out, stride, maxorig);
    } else
      throw std::runtime_error("invalid mode");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "spatial_dynamics: " << e.what() << '\n';
    return 1;
  }
}
