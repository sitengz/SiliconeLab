#ifndef PDMS_FILLER_COMPONENT_HPP
#define PDMS_FILLER_COMPONENT_HPP

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
#include <unordered_map>
#include <utility>
#include <vector>

namespace pdms_filler {

constexpr double kPi = 3.14159265358979323846;
constexpr double kBoundaryClearance = 1.0e-4;
constexpr double kDefaultMass = 74.0;
constexpr double kBondLength = 2.801;
constexpr double kBondAngleDegrees = 111.623;
constexpr double kEpsilon300K = 1.01287845;
constexpr double kSigma300K = 6.44584366;

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
    std::uint32_t seed = 20260727u;
    double minimum_separation = 4.5;
    double z_lower_fraction = -0.20;
    double z_upper_fraction = 0.20;
};

struct Atom {
    int molecule = 0;
    int type = 1;
    Vec3 position;
};

struct Bond {
    int type = 1;
    int a = 0;
    int b = 0;
};

struct Angle {
    int type = 1;
    int a = 0;
    int b = 0;
    int c = 0;
};

struct Dihedral {
    int type = 1;
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

inline double chain_mass(
    int length,
    double bead_mass = kDefaultMass
) {
    return length * bead_mass;
}

inline std::pair<Vec3, Vec3> perpendicular_basis(const Vec3& axis) {
    const Vec3 unit = normalized(axis);
    const Vec3 reference =
        std::fabs(unit.z) < 0.85 ? Vec3{0.0, 0.0, 1.0}
                                : Vec3{0.0, 1.0, 0.0};
    const Vec3 first = normalized(cross(reference, unit));
    return {first, normalized(cross(unit, first))};
}

inline bool candidate_allowed(
    const std::vector<Vec3>& placed,
    const Vec3& candidate
) {
    const int candidate_index = static_cast<int>(placed.size());
    for (int existing_index = 0; existing_index < candidate_index;
         ++existing_index) {
        const int path_length = candidate_index - existing_index;
        if (path_length <= 2) continue;
        const double minimum = path_length == 3 ? 4.0 : 4.5;
        if (norm2(candidate - placed[static_cast<std::size_t>(existing_index)]) <
            minimum * minimum)
            return false;
    }
    return true;
}

inline std::vector<Vec3> build_local_chain(
    int length,
    std::mt19937& random
) {
    std::vector<Vec3> chain(static_cast<std::size_t>(length));
    if (length == 1) return chain;
    chain[1] = {kBondLength, 0.0, 0.0};
    std::uniform_real_distribution<double> azimuth(0.0, 2.0 * kPi);

    for (int i = 2; i < length; ++i) {
        bool accepted = false;
        for (int attempt = 0; attempt < 200 && !accepted; ++attempt) {
            const Vec3 previous = normalized(
                chain[static_cast<std::size_t>(i - 1)] -
                chain[static_cast<std::size_t>(i - 2)]);
            const auto basis = perpendicular_basis(previous);
            const double alpha = kPi - radians(kBondAngleDegrees);
            const double phi = azimuth(random);
            const Vec3 direction =
                std::cos(alpha) * previous +
                std::sin(alpha) *
                    (std::cos(phi) * basis.first +
                     std::sin(phi) * basis.second);
            const Vec3 candidate =
                chain[static_cast<std::size_t>(i - 1)] +
                kBondLength * normalized(direction);
            std::vector<Vec3> placed(chain.begin(), chain.begin() + i);
            if (!candidate_allowed(placed, candidate)) continue;
            chain[static_cast<std::size_t>(i)] = candidate;
            accepted = true;
        }
        if (!accepted) return {};
    }

    Vec3 center;
    for (const Vec3& position : chain) center = center + position;
    center = center / static_cast<double>(length);
    for (Vec3& position : chain) position = position - center;
    return chain;
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
                                delta.z -=
                                    std::round(delta.z / box_.lz) * box_.lz;
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
        const auto coordinate = [](
            double value,
            double length,
            double width,
            int count
        ) {
            int index = static_cast<int>(
                std::floor((value + 0.5 * length) / width));
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

inline Component generate(
    const Settings& settings,
    const Box& box,
    const std::vector<Vec3>& existing_positions
) {
    if (settings.chains == 0) return {};
    if (settings.length <= 0)
        throw std::runtime_error("PDMS filler length must be positive");
    if (settings.minimum_separation <= 0.0)
        throw std::runtime_error("PDMS filler minimum separation must be positive");
    if (settings.z_lower_fraction < -0.5 ||
        settings.z_upper_fraction > 0.5 ||
        settings.z_lower_fraction >= settings.z_upper_fraction)
        throw std::runtime_error(
            "PDMS filler z-placement fractions must satisfy "
            "-0.5 <= lower < upper <= 0.5");

    Component component;
    component.atoms.reserve(
        static_cast<std::size_t>(settings.chains) *
        static_cast<std::size_t>(settings.length));
    CellList cells(box, settings.minimum_separation);
    for (const Vec3& position : existing_positions) cells.insert(position);
    std::mt19937 random(settings.seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (int molecule = 0; molecule < settings.chains; ++molecule) {
        std::vector<Vec3> accepted_positions;
        bool accepted = false;
        for (int attempt = 0; attempt < 1000 && !accepted; ++attempt) {
            const std::vector<Vec3> local =
                build_local_chain(settings.length, random);
            if (local.empty()) continue;
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
            for (const Vec3& position : local) {
                const Vec3 value = rotate(rotation, position);
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
                anchor_lower.x + unit(random) *
                    (anchor_upper.x - anchor_lower.x),
                anchor_lower.y + unit(random) *
                    (anchor_upper.y - anchor_lower.y),
                anchor_lower.z + unit(random) *
                    (anchor_upper.z - anchor_lower.z)
            };
            std::vector<Vec3> positions;
            positions.reserve(rotated.size());
            for (const Vec3& value : rotated)
                positions.push_back(anchor + value);
            if (cells.overlaps(positions)) continue;
            accepted_positions = std::move(positions);
            accepted = true;
        }
        if (!accepted) {
            std::ostringstream message;
            message << "Could not place PDMS filler chain " << molecule + 1
                    << " without overlap. Reduce filler loading, chain "
                       "length, minimum separation, or initial density.";
            throw std::runtime_error(message.str());
        }

        for (const Vec3& position : accepted_positions) cells.insert(position);
        const int first_atom = static_cast<int>(component.atoms.size()) + 1;
        for (const Vec3& position : accepted_positions)
            component.atoms.push_back({molecule + 1, 1, position});
        for (int i = 0; i + 1 < settings.length; ++i)
            component.bonds.push_back({1, first_atom + i, first_atom + i + 1});
        for (int i = 0; i + 2 < settings.length; ++i)
            component.angles.push_back(
                {1, first_atom + i, first_atom + i + 1, first_atom + i + 2});
        for (int i = 0; i + 3 < settings.length; ++i)
            component.dihedrals.push_back({
                1, first_atom + i, first_atom + i + 1,
                first_atom + i + 2, first_atom + i + 3
            });
    }
    return component;
}

inline PairParameters pair_parameters(double temperature) {
    if (std::fabs(temperature - 300.0) < 1.0e-12)
        return {kEpsilon300K, kSigma300K};
    const double shifted = temperature - 186.04682;
    const double denominator = 1.0 + std::exp(0.00758 * shifted);
    return {
        (4.77795 / denominator + 1.47169) * 0.350646,
        (7.86548e-05 * temperature + 1.27856) * 4.95013
    };
}

inline double repulsive_cutoff(const PairParameters& pair) {
    return pair.sigma * std::pow(2.0, 1.0 / 6.0);
}

inline double maximum_repulsive_cutoff(double temperature) {
    return repulsive_cutoff(pair_parameters(temperature));
}

inline void write_pair_matrix(
    std::ostream& output,
    double temperature,
    bool include_cutoff
) {
    const PairParameters pair = pair_parameters(temperature);
    output << std::fixed << std::setprecision(9);
    for (int i = 1; i <= 3; ++i) {
        for (int j = i; j <= 3; ++j) {
            output << "pair_coeff      " << i << ' ' << j << ' '
                   << pair.epsilon << ' ' << pair.sigma;
            if (include_cutoff)
                output << ' ' << repulsive_cutoff(pair);
            output << '\n';
        }
    }
}

} // namespace pdms_filler

#endif
