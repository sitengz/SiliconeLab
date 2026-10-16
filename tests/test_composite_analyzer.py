"""Composite snapshot contracts: scientific parity, NP isolation and validation."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

P = argparse.ArgumentParser()
P.add_argument("--generator", required=True, type=Path)
P.add_argument("--driver", required=True, type=Path)
P.add_argument("--backend-dir", required=True, type=Path)
ARGS, REST = P.parse_known_args()


class CompositeAnalyzerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def generate(self, base="v22", extra=""):
        native = ("strand_length=8\nstrand_count=4\ncrosslinker_length=8\nmoderator_count=0\n"
                  if base == "elastomer" else "n1=8\nm1=4\nn2=8\nm4=0\n")
        if "strand_topology=grafted" in extra:
            native = "strand_count=24\ncrosslinker_length=16\nmoderator_count=0\n"
        config = self.root / "case.conf"
        config.write_text(f"system=composite\nbase={base}\noutput_dir=case\n" + native +
                          "functionality=4\ndensity=0.1\ntarget_density=0.8\nnp_count=1\nnp_radius=10\nnp_bead_mass=148\n" + extra)
        p = subprocess.run([ARGS.generator, "--config", config], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return self.root / "case/data.case", self.root / "case/case.info"

    def analyze(self, data, info, name="analysis", extra=("--z1-mode", "export"), ok=True):
        out = self.root / name
        p = subprocess.run([sys.executable, ARGS.driver, data, info, "--backend-dir", ARGS.backend_dir,
                            "--output-dir", out, *extra], capture_output=True, text=True)
        if ok:
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr +
                             ((out / "analysis.log").read_text() if (out / "analysis.log").exists() else ""))
        return p, out

    def test_all_bases_matrix_parity_and_density(self):
        for base in ("v22", "v35", "elastomer"):
            with self.subTest(base=base):
                # Use a fresh directory for each native formulation.
                original = self.root
                self.root = original / base
                self.root.mkdir()
                data, info = self.generate(base)
                before = hashlib.sha256(data.read_bytes()).hexdigest()
                _, out = self.analyze(data, info)
                result = json.loads((out / "analysis.json").read_text())
                self.assertEqual(result["system"], "composite")
                quality = result["quality"]["directional_periodic_paths"]
                self.assertEqual(quality["status"], "provisional")
                self.assertFalse(quality["ordinary_connectivity_and_defect_counts_affected"])
                self.assertTrue(any("periodic self-image paths are provisional" in w for w in result["warnings"]))
                self.assertTrue(all(r["status"] == "complete" for r in result["analyses"][:10]))
                c = result["composite"]
                self.assertEqual(c["np_count"], 1)
                self.assertEqual(c["np_beads"], 13)
                self.assertAlmostEqual(c["np_mass_g_mol_equivalent"], 13 * 148)
                self.assertAlmostEqual(c["composite_density_g_cm3"],
                    c["matrix_mass_per_box_volume_g_cm3"] * c["total_mass_g_mol_equivalent"] / c["matrix_mass_g_mol_equivalent"])
                self.assertAlmostEqual(result["snapshot"]["density"], c["matrix_mass_per_box_volume_g_cm3"])
                self.assertEqual(before, hashlib.sha256(data.read_bytes()).hexdigest())
                precursor = self.root / "matrix.info"
                precursor.write_text(json.dumps(json.loads(info.read_text())["precursor_info"]))
                _, baseline = self.analyze(out / "matrix.data", precursor, name="baseline")
                for name in ("snapshot_statistics.tsv", "chemical_network_statistics.tsv",
                             "network_statistics.case.tsv", "strand_properties.case.tsv",
                             "effective_strands.case.tsv", "network_nodes.case.tsv",
                             "parent_molecules.case.tsv", "directional_self_paths.case.tsv"):
                    self.assertEqual((out / name).read_bytes(), (baseline / name).read_bytes(), name)
                self.root = original

    def test_generic_architectures(self):
        for topology in ("ring", "star", "grafted"):
            original = self.root
            self.root = original / topology
            self.root.mkdir()
            extra = f"strand_topology={topology}\n"
            if topology == "grafted":
                # Native graft generator requires length to be derived from its architecture.
                config_native = "backbone_length=12\nside_chain_length=4\ngraft_spacing=2\ngraft_functional_fraction=50\n"
                extra += config_native
            data, info = self.generate("elastomer", extra)
            _, out = self.analyze(data, info)
            self.assertEqual(json.loads((out / "analysis.json").read_text())["system"], "composite")
            self.root = original

    def test_oil_chemistry_and_legacy_format(self):
        data, info = self.generate("v35", "oil=copolymer\noil_length=4\nmps_percent=50\noil_wt=20\n")
        meta = json.loads(info.read_text())
        meta["format"] = "siliconelab-np-coating-info"
        info.write_text(json.dumps(meta))
        _, out = self.analyze(data, info)
        with (out / "oil_chain_properties.tsv").open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertTrue(rows)
        self.assertTrue(all(int(r["backbone_beads"]) == 4 for r in rows))
        self.assertEqual(sum(int(r["beads"]) for r in rows), meta["precursor_info"]["components"]["silicone_oil_filler"]["beads"])

    def test_composite_part2_part3_metadata_is_explicitly_unsupported(self):
        _, info = self.generate()
        spec = importlib.util.spec_from_file_location(
            "composite_analysis_common", ARGS.driver.parent / "siliconelab_analysis_common.py")
        common = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(common)
        meta = json.loads(info.read_text())
        for fmt in ("siliconelab-composite-info", "siliconelab-np-coating-info"):
            meta["format"] = fmt
            info.write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError, "Composite Part 2/Part 3 analysis is unsupported"):
                common.metadata(info)

    def test_known_reacted_graph_ignores_np_and_preserves_topology(self):
        data, info = self.generate("elastomer")
        meta = json.loads(info.read_text())["precursor_info"]
        strands, cross = meta["components"]["strands"], meta["components"]["crosslinkers"]
        ends = [[m * strands["N"] + s for s in meta["strand"]["reactive_bead_sites"]] for m in range(strands["M"])]
        sites = [[strands["beads"] + m * cross["N"] + s for s in meta["crosslinker"]["reactive_bead_sites"]]
                 for m in range(cross["M"])]
        pairs = [(ends[m][j], sites[j][m]) for m in range(strands["M"]) for j in (0, 1)]
        reacted = {a for pair in pairs for a in pair}
        lines = data.read_text().splitlines()
        bond_count = meta["topology"]["bonds"]
        section = ""
        for i, line in enumerate(lines):
            fields = line.split()
            if len(fields) == 2 and fields[1] == "bonds":
                lines[i] = f"{bond_count + len(pairs)} bonds"
            if line.startswith("Atoms"):
                section = "atoms"
                continue
            if line.startswith(("Bonds", "Angles", "Dihedrals", "Velocities")):
                section = "other"
            if section == "atoms" and fields and int(fields[0]) in reacted:
                fields[2] = "1"
                lines[i] = " ".join(fields)
        start = next(i for i, line in enumerate(lines) if line.startswith("Bonds")) + 1
        stop = next((i for i in range(start, len(lines)) if lines[i].startswith(("Angles", "Dihedrals", "Velocities"))), len(lines))
        lines[stop:stop] = [f"{bond_count+i+1} 2 {a} {b}" for i, (a, b) in enumerate(pairs)] + [""]
        data.write_text("\n".join(lines) + "\n")
        _, out = self.analyze(data, info)
        result = json.loads((out / "analysis.json").read_text())
        self.assertEqual(result["chemical_network"]["overall_conversion"], 1)
        self.assertEqual(result["network"]["active_strands"], 4)
        self.assertEqual(result["network"]["cycle_rank"], 3)
        self.assertGreater(result["network"]["affine_modulus_estimate"], 0)
        self.assertEqual(result["analyses"][10]["chains"], 4)
        mapping = json.loads((out / "z1/chain_mapping.json").read_text())
        matrix_atoms = json.loads(info.read_text())["precursor_beads"]
        self.assertTrue(all(int(a) <= matrix_atoms for r in mapping["chains"] for a in r["atom_ids"].split(",")))

    def test_np_bond_and_higher_topology_rejected(self):
        data, info = self.generate()
        original = data.read_text()
        np_id = json.loads(info.read_text())["precursor_beads"] + 1
        for i, (section, width) in enumerate((("Bonds", 4), ("Angles", 5), ("Dihedrals", 6))):
            lines = original.splitlines()
            start = lines.index(section) + 1
            row = next(j for j in range(start, len(lines)) if lines[j].strip())
            fields = lines[row].split()
            fields[-1] = str(np_id)
            lines[row] = " ".join(fields)
            data.write_text("\n".join(lines) + "\n")
            p, out = self.analyze(data, info, name=f"bad{i}", ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("NP covalent topology", p.stderr)
            self.assertFalse((out / "matrix.data").exists())

    def test_bad_metadata_np_atom_and_velocities_rejected(self):
        data, info = self.generate()
        original = info.read_text()
        for i, (key, value) in enumerate((("np_beads", 99), ("np_atom_type", 4),
                                          ("np_weight_percent", 99), ("mass_g_mol_per_bead", 74))):
            meta = json.loads(original)
            meta[key] = value
            info.write_text(json.dumps(meta))
            p, _ = self.analyze(data, info, name=f"bad{i}", ok=False)
            self.assertNotEqual(p.returncode, 0)
        info.write_text(original)
        with data.open("a") as f:
            f.write("\nVelocities\n\n1 0 0 0\n")
        p, _ = self.analyze(data, info, name="velocity", ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("incomplete composite velocities", p.stderr)

    def test_z1_np_obstacles_fail_closed_and_no_overwrite(self):
        data, info = self.generate()
        p, out = self.analyze(data, info, extra=(), ok=False)
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
        result = json.loads((out / "analysis.json").read_text())
        self.assertTrue(all(r["status"] == "complete" for r in result["analyses"][:10]))
        self.assertEqual(result["analyses"][11]["status"], "blocked")
        self.assertIn("NP obstacles", result["analyses"][11]["reason"])
        before = (out / "analysis.json").read_bytes()
        p, _ = self.analyze(data, info, ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertEqual(before, (out / "analysis.json").read_bytes())

    def test_np_velocities_do_not_change_matrix_temperature(self):
        data, info = self.generate("elastomer")
        meta = json.loads(info.read_text())
        original = data.read_text()
        temperatures = []
        for run, np_speed in enumerate((1, 10000)):
            rows = []
            for ident in range(1, meta["total_beads"] + 1):
                speed = np_speed if ident > meta["precursor_beads"] else (ident % 3 - 1) * 0.001
                rows.append(f"{ident} {speed} 0 0")
            data.write_text(original + "\nVelocities\n\n" + "\n".join(rows) + "\n")
            _, out = self.analyze(data, info, name=f"vel{run}")
            temperatures.append(json.loads((out / "analysis.json").read_text())["snapshot"]["velocity_temperature"])
        self.assertGreater(temperatures[0], 0)
        self.assertEqual(*temperatures)

    def test_np_atom_identity_and_nonfinite_values_rejected(self):
        data, info = self.generate()
        original = data.read_text().splitlines()
        start = next(i for i, line in enumerate(original) if line.startswith("Atoms")) + 1
        target = next(i for i in range(start, len(original)) if original[i].split() and
                      original[i].split()[0] == str(json.loads(info.read_text())["precursor_beads"] + 1))
        for index, (field, value) in enumerate(((1, "1"), (2, "1"), (3, "1"), (4, "nan"))):
            lines = list(original)
            row = lines[target].split()
            row[field] = value
            lines[target] = " ".join(row)
            data.write_text("\n".join(lines) + "\n")
            p, _ = self.analyze(data, info, name=f"identity{index}", ok=False)
            self.assertNotEqual(p.returncode, 0)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0], *REST])
