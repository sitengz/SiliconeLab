"""Validate a noncovalent rigid-NP composite and preserve its matrix subgraph.

The native analyzers consume the projected matrix with original atom/molecule
IDs, coordinates, velocities, box and covalent topology. No geometric repacking,
NP pseudo-bonds or changes to the scientific network reduction are introduced.
"""
import collections
import math
from pathlib import Path
import re

FORMATS = {"siliconelab-composite-info", "siliconelab-np-coating-info"}


def normalize_composite(info, normalize, geometry=None):
    if info.get("format_version") != 1 or info.get("representation") != "rigid_fcc":
        raise ValueError("composite analysis requires version-1 rigid_fcc metadata")
    if geometry not in (None, "bulk") or info.get("geometry") != "bulk":
        raise ValueError("composite structural analysis currently supports bulk only")
    matrix = normalize(info["precursor_info"])
    if matrix["source_system"] not in ("coating", "elastomer"):
        raise ValueError("composite matrix must be a supported elastomer or coating")
    if matrix["geometry"] != "bulk":
        raise ValueError("matrix and composite geometry disagree")
    keys = ("np_count", "beads_per_np", "np_beads", "precursor_beads", "total_beads",
            "precursor_molecules", "total_molecules", "np_molecule_id_start", "np_atom_type")
    for key in keys:
        if type(info.get(key)) is not int or info[key] <= 0:
            raise ValueError(f"invalid composite count: {key}")
    if (info["precursor_beads"] != matrix["composition"]["total_beads"] or
        info["precursor_molecules"] != matrix["composition"]["total_molecules"] or
        info["np_molecule_id_start"] != info["precursor_molecules"] + 1 or
        info["np_count"] * info["beads_per_np"] != info["np_beads"] or
        info["total_beads"] != info["precursor_beads"] + info["np_beads"] or
        info["total_molecules"] != info["precursor_molecules"] + info["np_count"]):
        raise ValueError("inconsistent composite/matrix counts or molecule ranges")
    if info["np_atom_type"] != max(map(int, matrix["expected_type_masses"])) + 1:
        raise ValueError("NP atom type must follow and be distinct from all matrix types")
    mass = info.get("mass_g_mol_per_bead")
    if not isinstance(mass, (int, float)) or not math.isfinite(mass) or mass <= 0:
        raise ValueError("invalid NP bead mass")
    if info.get("matrix_system", matrix["source_system"]) != matrix["source_system"]:
        raise ValueError("composite matrix system disagrees with nested precursor metadata")
    matrix.update(matrix_system=matrix["source_system"], matrix_source_format=matrix["source_format"],
                  source_system="composite", source_format=info["format"],
                  analysis_scope="matrix covalent network and matrix oil; NPs excluded",
                  composite={key: info[key] for key in keys})
    return matrix


def project_matrix(source, destination, info, normalized):
    """Fail closed on count/type/topology errors before writing a matrix data file."""
    sections = {}
    header = []
    section = None
    supported = {"Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"}
    coefficients = {"Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs",
                    "Dihedral Coeffs", "Improper Coeffs", "BondBond Coeffs", "BondAngle Coeffs"}
    with Path(source).open() as stream:
        title = next(stream, "")
        for line in stream:
            clean = line.split("#", 1)[0].strip()
            if not clean:
                continue
            if clean in supported | coefficients:
                if clean in sections:
                    raise ValueError("duplicate LAMMPS data section")
                section = clean
                sections[section] = []
            elif section is None:
                if any(x in clean for x in ("xy xz yz", "avec", "bvec", "cvec")):
                    raise ValueError("triclinic composite snapshots are unsupported")
                header.append(clean)
            else:
                if clean[0].isalpha():
                    raise ValueError(f"unsupported LAMMPS data section: {clean}")
                sections[section].append(clean.split())
    counts = {}
    bounds = {}
    for line in header:
        fields = line.split()
        if len(fields) == 2 and fields[1] in ("atoms", "bonds", "angles", "dihedrals", "impropers"):
            if fields[1] in counts:
                raise ValueError("duplicate snapshot count")
            counts[fields[1]] = int(fields[0])
        if len(fields) == 4 and fields[2] in ("xlo", "ylo", "zlo"):
            axis = fields[2][0]
            if axis in bounds or fields[3] != axis + "hi":
                raise ValueError("invalid box bounds")
            lo, hi = map(float, fields[:2])
            if not all(map(math.isfinite, (lo, hi))) or hi <= lo:
                raise ValueError("invalid box length")
            bounds[axis] = (lo, hi)
    if len(bounds) != 3 or counts.get("atoms") != info["total_beads"]:
        raise ValueError("composite atom count or box differs from metadata")
    type_header = [line for line in header if line.endswith(" atom types")]
    if len(type_header) != 1 or int(type_header[0].split()[0]) != info["np_atom_type"]:
        raise ValueError("composite atom-type count differs from metadata")
    masses = {}
    for fields in sections.get("Masses", []):
        if len(fields) != 2:
            raise ValueError("invalid mass row")
        typ, mass = int(fields[0]), float(fields[1])
        if typ in masses or not math.isfinite(mass) or mass <= 0:
            raise ValueError("duplicate or invalid mass")
        masses[typ] = mass
    expected = {int(k): float(v) for k, v in normalized["expected_type_masses"].items()}
    expected[info["np_atom_type"]] = info["mass_g_mol_per_bead"]
    if masses.keys() != expected.keys() or any(not math.isclose(masses[k], v, rel_tol=1e-8, abs_tol=1e-5)
                                             for k, v in expected.items()):
        raise ValueError("composite bead mass differs from metadata")
    atoms = {}
    np_ids = set()
    np_sizes = collections.Counter()
    matrix_mass = np_mass = 0.0
    np_masses = collections.defaultdict(float)
    for row in sections.get("Atoms", []):
        if len(row) not in (7, 10):
            raise ValueError("expected full-style atom row with optional image flags")
        ident, mol, typ = map(int, row[:3])
        values = list(map(float, row[3:7]))
        if not all(map(math.isfinite, values)):
            raise ValueError("nonfinite atom coordinates or charge")
        if len(row) == 10:
            list(map(int, row[7:]))
        if ident in atoms or not 1 <= ident <= info["total_beads"] or typ not in masses:
            raise ValueError("invalid or duplicate atom ID/type")
        atoms[ident] = row
        if typ == info["np_atom_type"]:
            if (ident <= info["precursor_beads"] or
                not info["np_molecule_id_start"] <= mol <= info["total_molecules"] or values[0] != 0):
                raise ValueError("NP atom identity, molecule range or charge mismatch")
            np_ids.add(ident)
            np_sizes[mol] += 1
            np_mass += masses[typ]
            np_masses[mol] += masses[typ]
        else:
            if ident > info["precursor_beads"] or not 1 <= mol <= info["precursor_molecules"]:
                raise ValueError("matrix atom identity or molecule range mismatch")
            matrix_mass += masses[typ]
    if (len(atoms) != info["total_beads"] or len(np_ids) != info["np_beads"] or
        len(np_sizes) != info["np_count"] or set(np_sizes.values()) != {info["beads_per_np"]}):
        raise ValueError("missing atoms or inconsistent per-NP bead counts")
    for label, field in (("matrix", "matrix_mass_g_mol_equivalent"), ("np", "np_mass_g_mol_equivalent"),
                         ("total", "total_mass_g_mol_equivalent")):
        actual = {"matrix": matrix_mass, "np": np_mass, "total": matrix_mass + np_mass}[label]
        if field in info and not math.isclose(actual, info[field], rel_tol=1e-7, abs_tol=1e-4):
            raise ValueError(f"{label} mass differs from composite metadata")
    weight = 100 * np_mass / (matrix_mass + np_mass)
    if not math.isclose(weight, info["np_weight_percent"], rel_tol=1e-7, abs_tol=1e-6):
        raise ValueError("NP weight percentage differs from snapshot")
    for name, width in (("Bonds", 4), ("Angles", 5), ("Dihedrals", 6), ("Impropers", 6)):
        rows = sections.get(name, [])
        if len(rows) != counts.get(name.lower(), 0):
            raise ValueError(f"{name} count differs from snapshot header")
        ids = set()
        for row in rows:
            if len(row) != width:
                raise ValueError(f"invalid {name} record")
            ident, typ, *ends = map(int, row)
            if ident <= 0 or ident in ids or typ <= 0 or any(a not in atoms for a in ends):
                raise ValueError(f"invalid or duplicate {name} record")
            ids.add(ident)
            if any(a in np_ids for a in ends):
                raise ValueError(f"NP covalent topology is unsupported: {name}")
    velocities = sections.get("Velocities", [])
    velocity_ids = set()
    for row in velocities:
        if len(row) != 4 or not all(math.isfinite(float(x)) for x in row[1:]):
            raise ValueError("invalid or nonfinite velocity")
        ident = int(row[0])
        if ident not in atoms or ident in velocity_ids:
            raise ValueError("invalid or duplicate velocity atom ID")
        velocity_ids.add(ident)
    if velocities and len(velocity_ids) != len(atoms):
        raise ValueError("incomplete composite velocities")
    # Force coefficients are unnecessary for snapshot graph/geometry analysis.
    # The original file remains the authoritative full simulation snapshot.
    with Path(destination).open("x") as stream:
        stream.write("SiliconeLab matrix projection; source IDs, box and covalent topology preserved\n\n")
        for line in header:
            if line.endswith(" atoms"):
                line = f'{info["precursor_beads"]} atoms'
            elif line.endswith(" atom types"):
                line = f'{info["np_atom_type"] - 1} atom types'
            stream.write(line + "\n")
        for name in ("Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers"):
            if name not in sections:
                continue
            stream.write("\n" + name + (" # full" if name == "Atoms" else "") + "\n\n")
            for row in sections[name]:
                ident = int(row[0])
                if name == "Masses" and ident == info["np_atom_type"]:
                    continue
                if name in ("Atoms", "Velocities") and ident in np_ids:
                    continue
                stream.write(" ".join(row) + "\n")
    volume = math.prod(hi - lo for lo, hi in bounds.values())
    return {
        "total_beads": len(atoms), "matrix_beads": len(atoms) - len(np_ids),
        "np_beads": len(np_ids), "np_count": len(np_sizes), "beads_per_np": info["beads_per_np"],
        "total_molecules": info["total_molecules"], "matrix_molecules": info["precursor_molecules"],
        "matrix_mass_g_mol_equivalent": matrix_mass, "np_mass_g_mol_equivalent": np_mass,
        "total_mass_g_mol_equivalent": matrix_mass + np_mass, "np_weight_percent": weight,
        "volume_A3": volume, "composite_density_g_cm3": (matrix_mass + np_mass) / (0.602214076 * volume),
        "matrix_mass_per_box_volume_g_cm3": matrix_mass / (0.602214076 * volume),
        "density_denominator": "full current composite box; no NP excluded-volume correction",
        "chemical_graph": "matrix only; rigid NPs have no covalent topology",
        "projection_coefficients": "omitted; static geometry and topology require no pair-force evaluation",
        "structural_modulus_scope": "matrix covalent-network estimate per full box volume; no reinforcement prediction",
        "np_bodies": [{"molecule": m, "beads": np_sizes[m], "mass_g_mol": np_masses[m]}
                      for m in sorted(np_sizes)],
    }
