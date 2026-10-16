import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("bulk_workflow", ROOT / "scripts/bulk_workflow.py")
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

NATIVE = """units real
boundary p p p
atom_style full
bond_style harmonic
angle_style hybrid harmonic quartic
dihedral_style nharmonic
special_bonds lj 0 0 0.5
comm_modify cutoff 60
pair_style lj/cut 8
read_data data.initial
bond_coeff 2 115 7.2
pair_coeff 1 1 .5 6.6 7.4
pair_style lj/gromacs 12 15
bond_coeff 2 115 2.8
pair_coeff 1 1 1.01 6.44
pair_coeff 1 6 .808 6.44
pair_coeff 6 6 .808 6.44
"""


class BulkContracts(unittest.TestCase):
    def test_actual_lammps_errors_are_distinct_from_warning_prose(self):
        self.assertFalse(w.has_lammps_error("WARNING: Communication cutoff may lead to errors.\n"))
        self.assertTrue(w.has_lammps_error("ERROR: Lost atoms\n"))
        self.assertTrue(w.has_lammps_error("  ERROR on proc 0: Missing bond atoms\n"))
        self.assertTrue(w.has_lammps_error("100 300 nan 1.2\n"))

    def test_seven_configs_and_fixed_matrix(self):
        with tempfile.TemporaryDirectory() as d:
            batch = Path(d) / "batch"
            plan = w.initialize(batch)
            self.assertEqual(len(plan["cases"]), 7)
            self.assertEqual([r["lambda"] for r in plan["cases"]], [None, 1, .9, .8, .7, .6, .5])
            for c in plan["cases"]:
                text = (batch / c["id"] / "configs/model.conf").read_text()
                for setting in ("m1 = 900", "n1 = 128", "n2 = 32", "functionality = 8", "m4 = 6", "tensile_erate = 2e-8", "tensile_strain = 1"):
                    self.assertIn(setting, text)
                if c["lambda"] is not None:
                    self.assertIn("np_count = 32", text)
                    self.assertIn("np_interaction = native_pdms", text)
            with self.assertRaises(FileExistsError):
                w.initialize(batch)

    def test_cold_header_uses_scaled_cold_pairs_and_bonds(self):
        text = w.cold_header(NATIVE, "../formation/data.final")
        self.assertIn("pair_style lj/gromacs 12 15", text)
        self.assertNotIn("lj/cut", text)
        self.assertIn("bond_coeff 2 115 2.8", text)
        self.assertNotIn("115 7.2", text)
        self.assertIn("pair_coeff 1 6 .808 6.44", text)

    def test_independent_xyz_protocol_and_rigid_ordering(self):
        for axis in "xyz":
            for composite in (False, True):
                t = w.tensile_input(NATIVE, "../formation/data.final", axis, composite, 2e-8, 1)
                self.assertIn("run 10000000", t)
                self.assertIn("read_data ../formation/data.final", t)
                self.assertIn(axis + " erate ${erate}", t)
                self.assertNotIn("velocity ", t)
                self.assertNotIn("bond/create", t)
                for transverse in set("xyz") - {axis}:
                    self.assertIn(transverse + " 1 1 500", t)
                if composite:
                    self.assertLess(t.index("fix np_integrate"), t.index("fix integrate precursor npt"))
                    self.assertIn("dilate precursor", t)
                self.assertIn("run 1 start 0 stop $(step+1:%.0f)", t)

    def test_submission_serializes_each_chain_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            batch = Path(d) / "batch"
            plan = w.initialize(batch)
            state = w.load(batch / "workflow-state.json")
            state["status"] = "validated"
            for c in plan["cases"]:
                for stage in w.STAGES:
                    state["stages"][c["id"] + ":" + stage] = {"directory": str(batch / c["id"] / stage), "job_id": None}
            w.save(batch / "workflow-state.json", state)
            submitted = []
            def output(args, **kwargs):
                if args[0] == "sbatch":
                    submitted.append(args)
                    return str(1000 + len(submitted)) + "\n"
                return ""
            with patch.object(w.subprocess, "check_output", side_effect=output), patch.dict(os.environ, USER="test"):
                w.submit(batch)
                w.submit(batch)
            self.assertEqual(len(submitted), 28)
            for i, args in enumerate(submitted):
                dependencies = [x for x in args if x.startswith("--dependency=")]
                self.assertEqual(dependencies, [] if i % 4 == 0 else ["--dependency=afterok:" + str(1000 + i)])
                if i % 4:
                    self.assertIn("--kill-on-invalid-dep=yes", args)

    def test_eight_independent_jobs_blocks_new_roots(self):
        with tempfile.TemporaryDirectory() as d:
            batch = Path(d) / "batch"
            plan = w.initialize(batch)
            state = w.load(batch / "workflow-state.json")
            state["status"] = "validated"
            first = plan["cases"][0]["id"] + ":formation"
            state["stages"][first] = {"directory": str(batch), "job_id": None}
            w.save(batch / "workflow-state.json", state)
            with patch.object(w.subprocess, "check_output", return_value="\n".join(f"{i}|(null)|other-{i}" for i in range(8))), patch.dict(os.environ, USER="test"):
                with self.assertRaisesRegex(RuntimeError, "Eight independent"):
                    w.submit(batch)


if __name__ == "__main__":
    unittest.main()
