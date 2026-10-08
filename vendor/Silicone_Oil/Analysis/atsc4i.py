#!/usr/bin/env python3
"""Calculate raw Mordred ATSC4i for one CG silicone-oil chain.

The first repeat's Si receives a methyl cap; the final repeat's O receives
a Si(CH3)3 cap: Me-[Si(R)(R')-O]n-SiMe3. CG backbone atom IDs increase from
the Si-side end to the O-side end in data files written by this generator.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


BACKBONE_TYPES = {1: "D", 4: "M"}


def data_path_from_config(config_path: Path) -> tuple[Path, dict[str, str]]:
    """Resolve the data file next to the config from an explicit output key."""
    settings: dict[str, str] = {}
    for line_number, raw in enumerate(config_path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"{config_path}:{line_number}: expected key = value")
        key, value = (part.strip() for part in line.split("=", 1))
        settings[key.replace("-", "_")] = value.strip("\"'")
    if "output" not in settings:
        raise ValueError(
            f"{config_path} has no output key; use --data with the generated data file"
        )
    requested = Path(settings["output"])
    case_name = requested.name.removeprefix("data.")
    if not case_name:
        raise ValueError(f"Cannot derive a case name from output = {requested}")
    return (requested if requested.is_absolute() else config_path.parent / requested), settings


def sequence_from_data(data_path: Path, molecule_id: int) -> str:
    """Extract one oriented D/M chain from this generator's LAMMPS data file."""
    if molecule_id < 1:
        raise ValueError("molecule ID must be positive")
    atom_types: dict[int, int] = {}
    backbone: dict[int, int] = {}
    adjacency: dict[int, set[int]] = defaultdict(set)
    pendant_neighbors: dict[int, int] = defaultdict(int)
    section = "header"
    saw_atoms = False
    saw_bonds = False

    with data_path.open() as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if line == "Atoms":
                section = "atoms"
                saw_atoms = True
                continue
            if line == "Bonds":
                section = "bonds"
                saw_bonds = True
                continue
            if section == "bonds" and line == "Angles":
                break
            if not line:
                continue
            fields = line.split()
            if section == "atoms":
                if len(fields) < 3:
                    raise ValueError(f"Malformed atom record in {data_path}: {line}")
                atom_id, mol_id, atom_type = map(int, fields[:3])
                if mol_id == molecule_id:
                    atom_types[atom_id] = atom_type
                    if atom_type in BACKBONE_TYPES:
                        backbone[atom_id] = atom_type
                    elif atom_type != 5:
                        raise ValueError(
                            f"Unsupported atom type {atom_type} in molecule {molecule_id}"
                        )
            elif section == "bonds":
                if len(fields) < 4:
                    raise ValueError(f"Malformed bond record in {data_path}: {line}")
                _, bond_type, left, right = map(int, fields[:4])
                left_in = left in atom_types
                right_in = right in atom_types
                if left_in != right_in:
                    raise ValueError("A bond crosses molecule boundaries")
                if not left_in:
                    continue
                if left in backbone and right in backbone:
                    if bond_type not in (1, 2):
                        raise ValueError("Backbone bond has an unexpected type")
                    adjacency[left].add(right)
                    adjacency[right].add(left)
                elif {atom_types[left], atom_types[right]} == {4, 5}:
                    if bond_type != 3:
                        raise ValueError("MPS pendant bond has an unexpected type")
                    mps_id = left if atom_types[left] == 4 else right
                    pendant_neighbors[mps_id] += 1
                else:
                    raise ValueError("Unexpected bond within selected molecule")

    if not saw_atoms or not saw_bonds:
        raise ValueError(f"Expected Atoms and Bonds sections in {data_path}")
    if not backbone:
        raise ValueError(f"Molecule {molecule_id} has no backbone in {data_path}")
    mps_count = sum(atom_type == 4 for atom_type in backbone.values())
    pendant_count = sum(atom_type == 5 for atom_type in atom_types.values())
    if pendant_count != mps_count or any(
        pendant_neighbors[atom_id] != 1
        for atom_id, atom_type in backbone.items() if atom_type == 4
    ):
        raise ValueError("MPS backbone/pendant counts or bonds do not match")

    if len(backbone) == 1:
        order = list(backbone)
    else:
        ends = [atom_id for atom_id in backbone if len(adjacency[atom_id]) == 1]
        if len(ends) != 2:
            raise ValueError("Selected backbone is not a simple linear chain")
        order = [min(ends)]
        previous = None
        while True:
            neighbors = adjacency[order[-1]] - ({previous} if previous else set())
            if not neighbors:
                break
            if len(neighbors) != 1 or len(order) >= len(backbone):
                raise ValueError("Selected backbone branches or contains a cycle")
            previous, next_atom = order[-1], next(iter(neighbors))
            order.append(next_atom)
    if len(order) != len(backbone) or order != sorted(order):
        raise ValueError(
            "Backbone IDs do not follow the generator's Si-to-O chain order"
        )
    return "".join(BACKBONE_TYPES[backbone[atom_id]] for atom_id in order)


def capped_smiles(sequence: str) -> str:
    """Build Me-[Si(Me)(R)-O]n-SiMe3 (R=Me for D, phenyl for M)."""
    sequence = sequence.upper()
    if not sequence or any(unit not in "DM" for unit in sequence):
        raise ValueError("Sequence must contain one or more D/M repeat symbols")
    repeat = {"D": "[Si](C)(C)O", "M": "[Si](C)(c1ccccc1)O"}
    return "C" + "".join(repeat[unit] for unit in sequence) + "[Si](C)(C)C"


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def calculate_atsc4i(sequence: str) -> dict[str, object]:
    try:
        from rdkit import Chem, __version__ as rdkit_version
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        from mordred import Calculator
        from mordred.Autocorrelation import ATSC
    except ImportError as error:
        raise RuntimeError(
            "RDKit and Mordred are required; install Analysis/requirements.txt"
        ) from error

    smiles = capped_smiles(sequence)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse reconstructed SMILES: {smiles}")
    with_hydrogens = Chem.AddHs(mol)
    descriptor = Calculator([ATSC(4, "i")], ignore_3D=True)(with_hydrogens)[0]
    try:
        value = float(descriptor)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"Mordred could not calculate ATSC4i: {descriptor}") from error
    if not math.isfinite(value):
        raise RuntimeError(f"Mordred returned non-finite ATSC4i: {descriptor}")
    return {
        "ATSC4i": value,
        "smiles": Chem.MolToSmiles(mol),
        "formula": CalcMolFormula(with_hydrogens),
        "atom_count_including_hydrogen": with_hydrogens.GetNumAtoms(),
        "rdkit_version": rdkit_version,
        "mordred_version": package_version("mordredcommunity") or package_version("mordred"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--config", type=Path, help="Generator config with output key")
    source.add_argument("--data", type=Path, help="Generated CG LAMMPS data file")
    source.add_argument("--sequence", help="D/M sequence ordered Si-side to O-side")
    parser.add_argument("--molecule-id", type=int, default=1)
    args = parser.parse_args()

    settings: dict[str, str] = {}
    data_path = args.data
    if args.config is not None:
        data_path, settings = data_path_from_config(args.config)
    if args.sequence is not None:
        sequence = args.sequence.upper()
    else:
        if not data_path.is_file():
            raise SystemExit(
                f"Missing generated CG data: {data_path}. Generate the case first."
            )
        sequence = sequence_from_data(data_path, args.molecule_id)
        if "length" in settings and len(sequence) != int(settings["length"]):
            raise SystemExit("Extracted chain length does not match the config")

    result = {
        "sequence": sequence,
        "chain_length": len(sequence),
        "dms_repeats": sequence.count("D"),
        "mps_repeats": sequence.count("M"),
        "terminal_model": "Me-[Si(R)(R')-O]n-SiMe3",
        "molecule_id": args.molecule_id if data_path is not None else None,
        "cg_data": str(data_path) if data_path is not None else None,
        **calculate_atsc4i(sequence),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        raise SystemExit(f"Error: {error}") from None
