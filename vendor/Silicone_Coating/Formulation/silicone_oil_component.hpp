#ifndef SILICONE_OIL_COMPONENT_HPP
#define SILICONE_OIL_COMPONENT_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace silicone_oil {

constexpr double kPi = 3.14159265358979323846;
constexpr double kBoundaryClearance = 1.0e-4;
constexpr double kDmsMass = 74.0;
constexpr double kMpsBackboneMass = 59.1204;
constexpr double kMpsPendantMass = 77.106;
constexpr double kMpsRepeatMass = kMpsBackboneMass + kMpsPendantMass;
constexpr double kDmsBondLength = 2.801;
constexpr double kMpsBackboneBondLength = 2.8039;
constexpr double kMpsPendantBondLength = 3.12497;
constexpr double kDmsAngleDegrees = 111.623;
constexpr double kMpsBackboneAngleDegrees = 110.566;
constexpr double kMpsPendantAngleDegrees = 110.746;
constexpr double kDmsEpsilon = 1.01287845;
constexpr double kDmsSigma = 6.44584366;
constexpr double kMpsBackboneEpsilon = 1.779864125;
constexpr double kMpsBackboneSigma = 5.475173;
constexpr double kMpsPendantEpsilon = 0.848858275;
constexpr double kMpsPendantSigma = 5.96750575;
constexpr double kMpsBackbonePendantEpsilon = 1.33135725;
constexpr double kMpsBackbonePendantSigma = 5.72729065;
constexpr double kDmsMpsEpsilonFactor = 0.579966;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3 operator*(double scale, const Vec3& value) {
    return {scale * value.x, scale * value.y, scale * value.z};
}

inline Vec3 operator/(const Vec3& value, double scale) {
    return {value.x / scale, value.y / scale, value.z / scale};
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

inline double norm2(const Vec3& value) {
    return dot(value, value);
}

inline double norm(const Vec3& value) {
    return std::sqrt(norm2(value));
}

inline Vec3 normalized(const Vec3& value) {
    const double length = norm(value);
    if (length < 1.0e-12)
        throw std::runtime_error("Cannot normalize a zero-length vector");
    return value / length;
}

inline double radians(double degrees) {
    return degrees * kPi / 180.0;
}

template <typename T>
inline T clamp_value(T value, T lower, T upper) {
    return std::max(lower, std::min(value, upper));
}

struct Box {
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    bool periodic_z = true;
};

struct Settings {
    int length = 0;
    int chains = 0;
    int mps_per_chain = 0;
    int chains_with_extra_mps = 0; // First this many chains use mps_per_chain+1.
    std::string sequence = "random";
    std::uint32_t seed = 20260727u;
    double minimum_separation = 4.5;
    double z_lower_fraction = 0.0;
    double z_upper_fraction = 0.5;
};

struct Atom {
    int molecule = 0;
    int type = 0;
    Vec3 position;
};

struct Bond {
    int type = 0;
    int a = 0;
    int b = 0;
};

struct Angle {
    int type = 0;
    int a = 0;
    int b = 0;
    int c = 0;
};

struct Dihedral {
    int type = 0;
    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
};

struct Component {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
};

struct PairParameters {
    double epsilon = 0.0;
    double sigma = 0.0;
};

struct LocalMolecule {
    std::vector<Vec3> backbone;
    std::vector<Vec3> pendants;
};

struct LocalSite {
    Vec3 position;
    int monomer = 0;
    bool pendant = false;
};

inline double chain_mass(
    int length,
    int mps_per_chain,
    double dms_mass = kDmsMass
) {
    return (length - mps_per_chain) * dms_mass +
           mps_per_chain * kMpsRepeatMass;
}

inline std::vector<bool> make_sequence(
    const Settings& settings,
    std::mt19937& random
) {
    std::vector<bool> is_mps(static_cast<std::size_t>(settings.length), false);
    if (settings.mps_per_chain == 0) return is_mps;
    if (settings.mps_per_chain == settings.length) {
        std::fill(is_mps.begin(), is_mps.end(), true);
        return is_mps;
    }
    if (settings.sequence == "random") {
        std::vector<int> sites(static_cast<std::size_t>(settings.length));
        for (int i = 0; i < settings.length; ++i)
            sites[static_cast<std::size_t>(i)] = i;
        std::shuffle(sites.begin(), sites.end(), random);
        for (int i = 0; i < settings.mps_per_chain; ++i)
            is_mps[static_cast<std::size_t>(sites[static_cast<std::size_t>(i)])] = true;
    } else if (settings.sequence == "alternating") {
        for (int i = 0; i < settings.length; ++i) {
            const int before = i * settings.mps_per_chain / settings.length;
            const int after = (i + 1) * settings.mps_per_chain / settings.length;
            if (after > before) is_mps[static_cast<std::size_t>(i)] = true;
        }
    } else if (settings.sequence == "block") {
        const int first = (settings.length - settings.mps_per_chain) / 2;
        for (int i = first; i < first + settings.mps_per_chain; ++i)
            is_mps[static_cast<std::size_t>(i)] = true;
    } else {
        throw std::runtime_error("Oil sequence must be random, alternating, or block");
    }
    return is_mps;
}

inline double backbone_bond_length(bool left_mps, bool right_mps) {
    return left_mps && right_mps ? kMpsBackboneBondLength : kDmsBondLength;
}

inline double backbone_angle_degrees(
    bool left_mps,
    bool center_mps,
    bool right_mps
) {
    return left_mps && center_mps && right_mps
        ? kMpsBackboneAngleDegrees
        : kDmsAngleDegrees;
}

inline std::pair<Vec3, Vec3> perpendicular_basis(const Vec3& axis) {
    const Vec3 unit = normalized(axis);
    const Vec3 reference =
        std::fabs(unit.z) < 0.85 ? Vec3{0.0, 0.0, 1.0}
                                : Vec3{0.0, 1.0, 0.0};
    const Vec3 first = normalized(cross(reference, unit));
    return {first, normalized(cross(unit, first))};
}

inline Vec3 random_unit_vector(std::mt19937& random) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double z = 2.0 * unit(random) - 1.0;
    const double phi = 2.0 * kPi * unit(random);
    const double radius = std::sqrt(std::max(0.0, 1.0 - z * z));
    return {radius * std::cos(phi), radius * std::sin(phi), z};
}

inline int topological_distance(const LocalSite& a, const LocalSite& b) {
    return std::abs(a.monomer - b.monomer) +
           static_cast<int>(a.pendant) + static_cast<int>(b.pendant);
}

inline bool local_pair_allowed(const LocalSite& a, const LocalSite& b) {
    const int path_length = topological_distance(a, b);
    if (path_length <= 2) return true;
    const double minimum = path_length == 3 ? 4.0 : 4.5;
    return norm2(a.position - b.position) >= minimum * minimum;
}

inline bool candidate_sites_allowed(
    const std::vector<LocalSite>& existing,
    const std::vector<LocalSite>& candidates
) {
    for (const LocalSite& candidate : candidates)
        for (const LocalSite& site : existing)
            if (!local_pair_allowed(candidate, site)) return false;
    for (std::size_t i = 0; i < candidates.size(); ++i)
        for (std::size_t j = i + 1; j < candidates.size(); ++j)
            if (!local_pair_allowed(candidates[i], candidates[j])) return false;
    return true;
}

inline Vec3 endpoint_pendant_position(
    const Vec3& center,
    const Vec3& neighbor,
    double phi
) {
    const Vec3 toward_neighbor = normalized(neighbor - center);
    const auto basis = perpendicular_basis(toward_neighbor);
    const double theta = radians(kMpsPendantAngleDegrees);
    const Vec3 direction =
        std::cos(theta) * toward_neighbor +
        std::sin(theta) *
            (std::cos(phi) * basis.first + std::sin(phi) * basis.second);
    return center + kMpsPendantBondLength * normalized(direction);
}

inline Vec3 interior_pendant_position(
    const Vec3& left,
    const Vec3& center,
    const Vec3& right,
    double sign
) {
    const Vec3 toward_left = normalized(left - center);
    const Vec3 toward_right = normalized(right - center);
    const Vec3 bisector = normalized(toward_left + toward_right);
    Vec3 normal = cross(toward_left, toward_right);
    normal = norm(normal) < 1.0e-10
        ? perpendicular_basis(bisector).first
        : normalized(normal);
    const double denominator = dot(bisector, toward_left);
    double scale = std::cos(radians(kMpsPendantAngleDegrees)) / denominator;
    scale = clamp_value(scale, -1.0, 1.0);
    const double normal_scale = std::sqrt(std::max(0.0, 1.0 - scale * scale));
    return center + kMpsPendantBondLength *
        normalized(scale * bisector + sign * normal_scale * normal);
}

inline LocalMolecule build_local_molecule(
    const std::vector<bool>& is_mps,
    std::mt19937& random
) {
    const int length = static_cast<int>(is_mps.size());
    LocalMolecule molecule;
    molecule.backbone.resize(static_cast<std::size_t>(length));
    molecule.pendants.resize(static_cast<std::size_t>(length));
    std::vector<LocalSite> placed;
    placed.reserve(static_cast<std::size_t>(
        length + std::count(is_mps.begin(), is_mps.end(), true)));
    std::uniform_real_distribution<double> azimuth(0.0, 2.0 * kPi);
    const int side_parity = std::uniform_int_distribution<int>(0, 1)(random);

    placed.push_back({molecule.backbone[0], 0, false});
    if (length == 1) {
        if (is_mps[0])
            molecule.pendants[0] =
                kMpsPendantBondLength * random_unit_vector(random);
    } else {
        molecule.backbone[1] = {
            backbone_bond_length(is_mps[0], is_mps[1]), 0.0, 0.0
        };
        placed.push_back({molecule.backbone[1], 1, false});
        if (is_mps[0]) {
            molecule.pendants[0] = endpoint_pendant_position(
                molecule.backbone[0], molecule.backbone[1], azimuth(random));
            placed.push_back({molecule.pendants[0], 0, true});
        }

        for (int i = 2; i < length; ++i) {
            bool accepted = false;
            for (int attempt = 0; attempt < 200 && !accepted; ++attempt) {
                const Vec3 previous = normalized(
                    molecule.backbone[static_cast<std::size_t>(i - 1)] -
                    molecule.backbone[static_cast<std::size_t>(i - 2)]);
                const auto basis = perpendicular_basis(previous);
                const double alpha = kPi - radians(backbone_angle_degrees(
                    is_mps[static_cast<std::size_t>(i - 2)],
                    is_mps[static_cast<std::size_t>(i - 1)],
                    is_mps[static_cast<std::size_t>(i)]));
                const double phi = azimuth(random);
                const Vec3 direction =
                    std::cos(alpha) * previous +
                    std::sin(alpha) *
                        (std::cos(phi) * basis.first +
                         std::sin(phi) * basis.second);
                const Vec3 backbone =
                    molecule.backbone[static_cast<std::size_t>(i - 1)] +
                    backbone_bond_length(
                        is_mps[static_cast<std::size_t>(i - 1)],
                        is_mps[static_cast<std::size_t>(i)]) *
                    normalized(direction);
                std::vector<LocalSite> candidates{{backbone, i, false}};
                Vec3 pendant;
                if (is_mps[static_cast<std::size_t>(i - 1)]) {
                    const double sign =
                        ((i - 1 + side_parity) % 2 == 0) ? 1.0 : -1.0;
                    pendant = interior_pendant_position(
                        molecule.backbone[static_cast<std::size_t>(i - 2)],
                        molecule.backbone[static_cast<std::size_t>(i - 1)],
                        backbone, sign);
                    candidates.push_back({pendant, i - 1, true});
                }
                if (!candidate_sites_allowed(placed, candidates)) continue;
                molecule.backbone[static_cast<std::size_t>(i)] = backbone;
                if (is_mps[static_cast<std::size_t>(i - 1)])
                    molecule.pendants[static_cast<std::size_t>(i - 1)] = pendant;
                placed.insert(placed.end(), candidates.begin(), candidates.end());
                accepted = true;
            }
            if (!accepted) return {};
        }

        if (is_mps[static_cast<std::size_t>(length - 1)]) {
            bool accepted = false;
            for (int attempt = 0; attempt < 200 && !accepted; ++attempt) {
                const Vec3 pendant = endpoint_pendant_position(
                    molecule.backbone[static_cast<std::size_t>(length - 1)],
                    molecule.backbone[static_cast<std::size_t>(length - 2)],
                    azimuth(random));
                const std::vector<LocalSite> candidates{
                    {pendant, length - 1, true}
                };
                if (!candidate_sites_allowed(placed, candidates)) continue;
                molecule.pendants[static_cast<std::size_t>(length - 1)] = pendant;
                placed.push_back(candidates.front());
                accepted = true;
            }
            if (!accepted) return {};
        }
    }

    Vec3 center;
    int atom_count = 0;
    for (int i = 0; i < length; ++i) {
        center = center + molecule.backbone[static_cast<std::size_t>(i)];
        ++atom_count;
        if (is_mps[static_cast<std::size_t>(i)]) {
            center = center + molecule.pendants[static_cast<std::size_t>(i)];
            ++atom_count;
        }
    }
    center = center / atom_count;
    for (int i = 0; i < length; ++i) {
        molecule.backbone[static_cast<std::size_t>(i)] =
            molecule.backbone[static_cast<std::size_t>(i)] - center;
        if (is_mps[static_cast<std::size_t>(i)])
            molecule.pendants[static_cast<std::size_t>(i)] =
                molecule.pendants[static_cast<std::size_t>(i)] - center;
    }
    return molecule;
}

inline std::vector<LocalSite> local_sites(
    const LocalMolecule& molecule,
    const std::vector<bool>& is_mps
) {
    std::vector<LocalSite> sites;
    for (std::size_t i = 0; i < molecule.backbone.size(); ++i)
        sites.push_back({molecule.backbone[i], static_cast<int>(i), false});
    for (std::size_t i = 0; i < molecule.pendants.size(); ++i)
        if (is_mps[i])
            sites.push_back({molecule.pendants[i], static_cast<int>(i), true});
    return sites;
}

inline std::array<double, 9> random_rotation(std::mt19937& random) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double u1 = unit(random);
    const double u2 = unit(random);
    const double u3 = unit(random);
    const double qx = std::sqrt(1.0 - u1) * std::sin(2.0 * kPi * u2);
    const double qy = std::sqrt(1.0 - u1) * std::cos(2.0 * kPi * u2);
    const double qz = std::sqrt(u1) * std::sin(2.0 * kPi * u3);
    const double qw = std::sqrt(u1) * std::cos(2.0 * kPi * u3);
    return {
        1.0 - 2.0 * (qy*qy + qz*qz), 2.0 * (qx*qy - qz*qw),
        2.0 * (qx*qz + qy*qw), 2.0 * (qx*qy + qz*qw),
        1.0 - 2.0 * (qx*qx + qz*qz), 2.0 * (qy*qz - qx*qw),
        2.0 * (qx*qz - qy*qw), 2.0 * (qy*qz + qx*qw),
        1.0 - 2.0 * (qx*qx + qy*qy)
    };
}

inline Vec3 rotate(const std::array<double, 9>& matrix, const Vec3& value) {
    return {
        matrix[0]*value.x + matrix[1]*value.y + matrix[2]*value.z,
        matrix[3]*value.x + matrix[4]*value.y + matrix[5]*value.z,
        matrix[6]*value.x + matrix[7]*value.y + matrix[8]*value.z
    };
}

class CellList {
public:
    CellList(const Box& box, double cutoff)
        : box_(box), cutoff2_(cutoff * cutoff) {
        nx_ = std::max(1, static_cast<int>(std::floor(box.lx / cutoff)));
        ny_ = std::max(1, static_cast<int>(std::floor(box.ly / cutoff)));
        nz_ = std::max(1, static_cast<int>(std::floor(box.lz / cutoff)));
        wx_ = box.lx / nx_;
        wy_ = box.ly / ny_;
        wz_ = box.lz / nz_;
    }

    void insert(const Vec3& position) {
        const auto cell = coordinates(position);
        cells_[key(cell[0], cell[1], cell[2])].push_back(position);
    }

    bool overlaps(const std::vector<Vec3>& positions) const {
        for (const Vec3& position : positions) {
            const auto cell = coordinates(position);
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        const int cx = wrapped(cell[0] + dx, nx_);
                        const int cy = wrapped(cell[1] + dy, ny_);
                        int cz = cell[2] + dz;
                        if (box_.periodic_z) {
                            cz = wrapped(cz, nz_);
                        } else if (cz < 0 || cz >= nz_) {
                            continue;
                        }
                        const auto found = cells_.find(key(cx, cy, cz));
                        if (found == cells_.end()) continue;
                        for (const Vec3& existing : found->second) {
                            Vec3 delta = position - existing;
                            delta.x -= std::round(delta.x / box_.lx) * box_.lx;
                            delta.y -= std::round(delta.y / box_.ly) * box_.ly;
                            if (box_.periodic_z)
                                delta.z -= std::round(delta.z / box_.lz) * box_.lz;
                            if (norm2(delta) < cutoff2_) return true;
                        }
                    }
                }
            }
        }
        return false;
    }

private:
    static int wrapped(int index, int count) {
        index %= count;
        return index < 0 ? index + count : index;
    }

    std::array<int, 3> coordinates(const Vec3& position) const {
        const auto coordinate = [](double value, double length, double width, int count) {
            int index = static_cast<int>(std::floor((value + 0.5 * length) / width));
            return clamp_value(index, 0, count - 1);
        };
        return {
            coordinate(position.x, box_.lx, wx_, nx_),
            coordinate(position.y, box_.ly, wy_, ny_),
            coordinate(position.z, box_.lz, wz_, nz_)
        };
    }

    long long key(int x, int y, int z) const {
        return (static_cast<long long>(x) * ny_ + y) * nz_ + z;
    }

    Box box_;
    double cutoff2_;
    int nx_ = 1;
    int ny_ = 1;
    int nz_ = 1;
    double wx_ = 0.0;
    double wy_ = 0.0;
    double wz_ = 0.0;
    std::unordered_map<long long, std::vector<Vec3>> cells_;
};

inline void add_topology(
    Component& component,
    const std::vector<bool>& is_mps,
    const std::vector<int>& backbone_ids,
    const std::vector<int>& pendant_ids
) {
    const int length = static_cast<int>(is_mps.size());
    for (int i = 0; i + 1 < length; ++i) {
        // Bond type 2 is reserved exclusively for formulation crosslinks.
        const int type = is_mps[static_cast<std::size_t>(i)] &&
                         is_mps[static_cast<std::size_t>(i + 1)] ? 3 : 1;
        component.bonds.push_back({
            type, backbone_ids[static_cast<std::size_t>(i)],
            backbone_ids[static_cast<std::size_t>(i + 1)]
        });
    }
    for (int i = 0; i < length; ++i) {
        if (!is_mps[static_cast<std::size_t>(i)]) continue;
        component.bonds.push_back({
            4, backbone_ids[static_cast<std::size_t>(i)],
            pendant_ids[static_cast<std::size_t>(i)]
        });
    }

    for (int i = 1; i + 1 < length; ++i) {
        const int type =
            is_mps[static_cast<std::size_t>(i - 1)] &&
            is_mps[static_cast<std::size_t>(i)] &&
            is_mps[static_cast<std::size_t>(i + 1)] ? 2 : 1;
        component.angles.push_back({
            type, backbone_ids[static_cast<std::size_t>(i - 1)],
            backbone_ids[static_cast<std::size_t>(i)],
            backbone_ids[static_cast<std::size_t>(i + 1)]
        });
    }
    for (int i = 0; i < length; ++i) {
        if (!is_mps[static_cast<std::size_t>(i)]) continue;
        if (i > 0)
            component.angles.push_back({
                3, backbone_ids[static_cast<std::size_t>(i - 1)],
                backbone_ids[static_cast<std::size_t>(i)],
                pendant_ids[static_cast<std::size_t>(i)]
            });
        if (i + 1 < length)
            component.angles.push_back({
                3, pendant_ids[static_cast<std::size_t>(i)],
                backbone_ids[static_cast<std::size_t>(i)],
                backbone_ids[static_cast<std::size_t>(i + 1)]
            });
    }

    for (int i = 0; i + 1 < length; ++i) {
        struct EndAtom { int id; bool pendant; bool mps_backbone; };
        std::vector<EndAtom> left;
        std::vector<EndAtom> right;
        if (i > 0)
            left.push_back({
                backbone_ids[static_cast<std::size_t>(i - 1)], false,
                is_mps[static_cast<std::size_t>(i - 1)]
            });
        if (is_mps[static_cast<std::size_t>(i)])
            left.push_back({
                pendant_ids[static_cast<std::size_t>(i)], true, false
            });
        if (i + 2 < length)
            right.push_back({
                backbone_ids[static_cast<std::size_t>(i + 2)], false,
                is_mps[static_cast<std::size_t>(i + 2)]
            });
        if (is_mps[static_cast<std::size_t>(i + 1)])
            right.push_back({
                pendant_ids[static_cast<std::size_t>(i + 1)], true, false
            });

        for (const EndAtom& a : left) {
            for (const EndAtom& d : right) {
                const int pendant_count =
                    static_cast<int>(a.pendant) + static_cast<int>(d.pendant);
                int type = 1;
                if (pendant_count == 2) type = 4;
                else if (pendant_count == 1) type = 3;
                else if (a.mps_backbone &&
                         is_mps[static_cast<std::size_t>(i)] &&
                         is_mps[static_cast<std::size_t>(i + 1)] &&
                         d.mps_backbone) type = 2;
                component.dihedrals.push_back({
                    type, a.id,
                    backbone_ids[static_cast<std::size_t>(i)],
                    backbone_ids[static_cast<std::size_t>(i + 1)], d.id
                });
            }
        }
    }
}

inline Component generate(
    const Settings& settings,
    const Box& box,
    const std::vector<Vec3>& existing_positions
) {
    if (settings.chains == 0) return {};
    if (settings.length <= 0 ||
        settings.mps_per_chain < 0 ||
        settings.mps_per_chain > settings.length)
        throw std::runtime_error("Invalid silicone-oil chain composition");
    if (settings.chains_with_extra_mps < 0 ||
        settings.chains_with_extra_mps >= settings.chains ||
        (settings.chains_with_extra_mps > 0 &&
         settings.mps_per_chain >= settings.length))
        throw std::runtime_error("Invalid silicone-oil chain composition distribution");
    if (settings.z_lower_fraction < -0.5 ||
        settings.z_upper_fraction > 0.5 ||
        settings.z_lower_fraction >= settings.z_upper_fraction)
        throw std::runtime_error(
            "Silicone-oil z-placement fractions must satisfy "
            "-0.5 <= lower < upper <= 0.5");

    Component component;
    component.atoms.reserve(static_cast<std::size_t>(settings.chains) *
        static_cast<std::size_t>(settings.length + settings.mps_per_chain) +
        static_cast<std::size_t>(settings.chains_with_extra_mps));
    CellList cells(box, settings.minimum_separation);
    for (const Vec3& position : existing_positions) cells.insert(position);
    std::mt19937 random(settings.seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (int molecule = 0; molecule < settings.chains; ++molecule) {
        Settings chain_settings = settings;
        if (molecule < settings.chains_with_extra_mps)
            ++chain_settings.mps_per_chain;
        const std::vector<bool> is_mps = make_sequence(chain_settings, random);
        std::vector<Vec3> accepted_positions;
        bool accepted = false;
        for (int attempt = 0; attempt < 1000 && !accepted; ++attempt) {
            const LocalMolecule local = build_local_molecule(is_mps, random);
            if (local.backbone.empty()) continue;
            const std::vector<LocalSite> sites = local_sites(local, is_mps);
            const auto rotation = random_rotation(random);
            std::vector<Vec3> rotated;
            Vec3 lower{
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()
            };
            Vec3 upper{
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest()
            };
            for (const LocalSite& site : sites) {
                const Vec3 value = rotate(rotation, site.position);
                rotated.push_back(value);
                lower.x = std::min(lower.x, value.x);
                lower.y = std::min(lower.y, value.y);
                lower.z = std::min(lower.z, value.z);
                upper.x = std::max(upper.x, value.x);
                upper.y = std::max(upper.y, value.y);
                upper.z = std::max(upper.z, value.z);
            }

            const Vec3 anchor_lower{
                -0.5*box.lx + kBoundaryClearance - lower.x,
                -0.5*box.ly + kBoundaryClearance - lower.y,
                settings.z_lower_fraction*box.lz +
                    kBoundaryClearance - lower.z
            };
            const Vec3 anchor_upper{
                0.5*box.lx - kBoundaryClearance - upper.x,
                0.5*box.ly - kBoundaryClearance - upper.y,
                settings.z_upper_fraction*box.lz -
                    kBoundaryClearance - upper.z
            };
            if (anchor_lower.x > anchor_upper.x ||
                anchor_lower.y > anchor_upper.y ||
                anchor_lower.z > anchor_upper.z)
                continue;
            const Vec3 anchor{
                anchor_lower.x + unit(random) * (anchor_upper.x - anchor_lower.x),
                anchor_lower.y + unit(random) * (anchor_upper.y - anchor_lower.y),
                anchor_lower.z + unit(random) * (anchor_upper.z - anchor_lower.z)
            };
            std::vector<Vec3> positions;
            for (const Vec3& value : rotated)
                positions.push_back(anchor + value);
            if (cells.overlaps(positions)) continue;
            accepted_positions = std::move(positions);
            accepted = true;
        }
        if (!accepted) {
            std::ostringstream message;
            message << "Could not place oil chain " << molecule + 1
                    << " without overlap. Reduce oil loading, chain length, "
                       "minimum separation, or initial density.";
            throw std::runtime_error(message.str());
        }

        for (const Vec3& position : accepted_positions) cells.insert(position);
        std::vector<int> backbone_ids(static_cast<std::size_t>(settings.length));
        std::vector<int> pendant_ids(static_cast<std::size_t>(settings.length), 0);
        std::size_t site = 0;
        for (int i = 0; i < settings.length; ++i, ++site) {
            const int id = static_cast<int>(component.atoms.size()) + 1;
            backbone_ids[static_cast<std::size_t>(i)] = id;
            component.atoms.push_back({
                molecule + 1,
                is_mps[static_cast<std::size_t>(i)] ? 4 : 1,
                accepted_positions[site]
            });
        }
        for (int i = 0; i < settings.length; ++i) {
            if (!is_mps[static_cast<std::size_t>(i)]) continue;
            const int id = static_cast<int>(component.atoms.size()) + 1;
            pendant_ids[static_cast<std::size_t>(i)] = id;
            component.atoms.push_back({
                molecule + 1, 5, accepted_positions[site++]
            });
        }
        add_topology(component, is_mps, backbone_ids, pendant_ids);
    }
    return component;
}

inline PairParameters dms_pair_parameters(double temperature) {
    if (std::fabs(temperature - 300.0) < 1.0e-12)
        return {kDmsEpsilon, kDmsSigma};
    const double shifted = temperature - 186.04682;
    const double denominator = 1.0 + std::exp(0.00758 * shifted);
    return {
        (4.77795 / denominator + 1.47169) * 0.350646,
        (7.86548e-05 * temperature + 1.27856) * 4.95013
    };
}

inline PairParameters pair_parameters(
    int type_i,
    int type_j,
    double temperature
) {
    if (type_i > type_j) std::swap(type_i, type_j);
    const PairParameters dms = dms_pair_parameters(temperature);
    const double epsilon_scale = dms.epsilon / kDmsEpsilon;
    const double sigma_scale = dms.sigma / kDmsSigma;
    const PairParameters backbone{
        kMpsBackboneEpsilon * epsilon_scale,
        kMpsBackboneSigma * sigma_scale
    };
    const PairParameters pendant{
        kMpsPendantEpsilon * epsilon_scale,
        kMpsPendantSigma * sigma_scale
    };
    const PairParameters backbone_pendant{
        kMpsBackbonePendantEpsilon * epsilon_scale,
        kMpsBackbonePendantSigma * sigma_scale
    };
    if (type_i <= 3 && type_j <= 3) return dms;
    if (type_i <= 3 && type_j == 4)
        return {kDmsMpsEpsilonFactor *
                    std::sqrt(dms.epsilon * backbone.epsilon),
                0.5 * (dms.sigma + backbone.sigma)};
    if (type_i <= 3 && type_j == 5)
        return {kDmsMpsEpsilonFactor *
                    std::sqrt(dms.epsilon * pendant.epsilon),
                0.5 * (dms.sigma + pendant.sigma)};
    if (type_i == 4 && type_j == 4) return backbone;
    if (type_i == 4 && type_j == 5) return backbone_pendant;
    if (type_i == 5 && type_j == 5) return pendant;
    throw std::logic_error("Unsupported silicone-oil atom-type pair");
}

inline double repulsive_cutoff(const PairParameters& pair) {
    return pair.sigma * std::pow(2.0, 1.0 / 6.0);
}

inline double maximum_repulsive_cutoff(double temperature) {
    double maximum = 0.0;
    for (int i = 1; i <= 5; ++i)
        for (int j = i; j <= 5; ++j)
            maximum = std::max(
                maximum, repulsive_cutoff(pair_parameters(i, j, temperature)));
    return maximum;
}

inline void write_pair_matrix(
    std::ostream& output,
    double temperature,
    bool include_cutoff
) {
    output << std::fixed << std::setprecision(9);
    for (int i = 1; i <= 5; ++i) {
        for (int j = i; j <= 5; ++j) {
            const PairParameters pair = pair_parameters(i, j, temperature);
            output << "pair_coeff      " << i << ' ' << j << ' '
                   << pair.epsilon << ' ' << pair.sigma;
            if (include_cutoff)
                output << ' ' << repulsive_cutoff(pair);
            output << '\n';
        }
    }
}

} // namespace silicone_oil

#endif
