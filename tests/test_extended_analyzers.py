#!/usr/bin/env python3
"""Independent analytic fixtures for surface mechanics, dynamics and missing files."""

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import siliconelab_part3 as p3
from siliconelab_analysis_common import metadata, write_json, read_tsv

ARGS = None


class ExtendedTests(unittest.TestCase):
    def test_surface_units_blocks_and_wall_invalidity(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            p = out / "pressure.dat"
            # 100 A cell, 2 atm normal-minus-lateral -> 1.01325 mN/m.
            p.write_text(
                "".join(
                    f"{i*1000000} 300 -100 1 1 3 50 50 100 0 0\n" for i in range(1, 51)
                )
            )
            _, blocks, r = p3.surface_table(p, out)
            self.assertEqual(r["completed_blocks"], 10)
            self.assertEqual([b["samples"] for b in blocks], [5] * 10)
            self.assertAlmostEqual(r["block_mean_mN_m"], 1.01325)
            self.assertEqual(r["status"], "complete")
            p.write_text(p.read_text().replace("100 0 0", "100 1 0"))
            _, _, r = p3.surface_table(p, out)
            self.assertEqual(r["status"], "invalid_wall_contact")

    def test_tensile_vacuum_correction(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            p = out / "stress.dat"
            lines = []
            for i in range(1, 201):
                e = i * 0.001
                difference = 2 * e
                px = -difference / (0.101325 * 2)
                lines.append(
                    f"{i*1000} {i*5000} 300 {e} 0 0 {px} 0 0 50 50 100 250000\n"
                )
            p.write_text("".join(lines))
            _, r = p3.tensile_table(p, "x", 50, out)
            self.assertEqual(r["status"], "complete")
            self.assertAlmostEqual(
                r["small_strain_stress_difference_slope_MPa"]["slope"], 2
            )

    def test_phase_reset_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "p"
            p.write_text(
                "5000 300 0 1 1 3 50 50 100 0 0\n5000 300 0 1 1 3 50 50 100 0 0\n"
            )
            with self.assertRaises(ValueError):
                p3.surface_table(p, Path(d))

    def test_streaming_unwrap_drift_and_nonuniform_lags(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            config = root / "oil.conf"
            case = root / "tiny"
            config.write_text(
                f"system = oil\noutput_dir = {case}\nchain_count = 4 2\nmps_percent = 0\n"
            )
            subprocess.run(
                [ARGS.generator, "--config", str(config)],
                check=True,
                capture_output=True,
            )
            data = case / "data.tiny"
            source, norm = metadata(case / "tiny.info", "film")
            write_json(root / "normalized.info", norm)
            # Extract atom identities and manufacture a small, exactly solvable trajectory.
            atoms = []
            inside = False
            for line in data.read_text().splitlines():
                w = line.split()
                if w and w[0] == "Atoms":
                    inside = True
                    continue
                if inside and w:
                    if not w[0].isdigit():
                        break
                    atoms.append((int(w[0]), int(w[1]), int(w[2])))
            dump = root / "dump"
            text = ""
            for step in [0, 1000, 3000]:
                text += f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n8\nITEM: BOX BOUNDS pp pp ff\n0 10\n0 10\n-50 50\nITEM: ATOMS id mol type x y z ix iy iz\n"
                for atom, mol, typ in atoms:
                    # Opposite motions plus a common translation; crossing x must unwrap.
                    xu = (step / 1000) * (1 if mol == 1 else -1) + step / 1000 * 2
                    x = xu % 10
                    ix = int((xu - x) / 10)
                    text += f"{atom} {mol} {typ} {x} 0 0 {ix} 0 0\n"
            dump.write_text(text)
            out = root / "results"
            subprocess.run(
                [
                    ARGS.backend,
                    "dynamics",
                    str(data),
                    str(root / "normalized.info"),
                    str(out),
                    str(dump),
                    "1",
                    "3",
                ],
                check=True,
                capture_output=True,
            )
            rows = read_tsv(out / "dynamics.tsv")
            r = next(
                r
                for r in rows
                if r["lag_steps"] == "3000"
                and r["component"] == "oil"
                and r["particle"] == "molecular_COM"
                and r["layer_center_relative_A"] == "global"
            )
            self.assertAlmostEqual(float(r["MSD_x_A2"]), 9)
            self.assertEqual(
                json.loads((out / "sampling.json").read_text())["frames"], 3
            )
            dump.write_text(text.replace("3000\nITEM: NUMBER", "1000\nITEM: NUMBER"))
            p = subprocess.run(
                [
                    ARGS.backend,
                    "dynamics",
                    str(data),
                    str(root / "normalized.info"),
                    str(out),
                    str(dump),
                    "1",
                    "3",
                ],
                capture_output=True,
            )
            self.assertNotEqual(p.returncode, 0)

    def test_part2_four_schema_outputs(self):
        if importlib.util.find_spec("numpy") is None:
            self.skipTest("NumPy unavailable in this runtime")
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for name in ["oil", "elastomer", "coating-v22", "coating-v35"]:
                src = (ROOT / "examples" / f"{name}.conf").read_text()
                lines = [
                    line
                    for line in src.splitlines()
                    if not line.startswith("output_dir")
                ]
                cfg = root / (name + ".conf")
                case = root / name
                cfg.write_text("\n".join(lines) + "\noutput_dir = " + str(case) + "\n")
                subprocess.run(
                    [ARGS.generator, "--config", str(cfg)],
                    check=True,
                    capture_output=True,
                )
                out = root / (name + "-analysis")
                subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / "src/siliconelab_part2.py"),
                        str(case / ("data." + name)),
                        str(case / (name + ".info")),
                        "--output-dir",
                        str(out),
                        "--backend-dir",
                        str(Path(ARGS.backend).parent),
                    ],
                    check=True,
                    capture_output=True,
                )
                r = json.loads((out / "analysis.json").read_text())
                self.assertEqual(len(r["analyses"]), 12)
                self.assertEqual(r["status"], "complete")
                self.assertEqual(r["analyses"][7]["status"], "missing_output")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--generator", default=str(ROOT / "build/siliconelab_generator"))
    p.add_argument(
        "--backend", default=str(ROOT / "build/siliconelab_spatial_dynamics")
    )
    ARGS, rest = p.parse_known_args()
    unittest.main(argv=[sys.argv[0]] + rest)
