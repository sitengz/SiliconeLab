"""Interface geometry, graph, provenance and continuation contracts."""

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import siliconelab_interface as w

spec = importlib.util.spec_from_file_location(
    "composite_dynamics", ROOT / "scripts/composite_dynamics.py"
)
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
p = argparse.ArgumentParser()
p.add_argument("--generator", type=Path, required=True)
p.add_argument("--backend-dir", type=Path, required=True)
ARGS, REST = p.parse_known_args()


class InterfaceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def configuration(self, **extra):
        p = self.root / "analysis.conf"
        p.write_text(
            json.dumps(
                dict(data="source.data", info="source.info", output_dir="out", **extra)
            )
        )
        return w.configure(p)

    def test_config_rejects_unknown_invalid_and_unresolved_boundaries(self):
        c = self.configuration()
        self.assertEqual(c["interface_gap_A"], 10)
        self.assertTrue(Path(c["data"]).is_absolute())
        for extra in (
            {"lambda": 1},
            {"volume_points": True},
            {"interface_gap_A": 12},
            {"reference_gap_A": 5},
            {"shell_edges_A": [0, 10, 5]},
            {"max_origins": 0},
        ):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                self.configuration(**extra)

    def test_parallel_edges_survive_but_dangling_tree_is_pruned(self):
        edges = [
            dict(id=i + 1, u=a, v=b, status="active")
            for i, (a, b) in enumerate([(1, 2), (2, 3), (3, 1), (3, 4)])
        ]
        core, giant = w.pruned_core(edges)
        self.assertEqual(core, {1, 2, 3})
        self.assertEqual(giant, {1, 2, 3, 4})
        core, _ = w.pruned_core(
            [
                dict(id=1, u=1, v=2, status="active"),
                dict(id=2, u=1, v=2, status="active"),
            ]
        )
        self.assertEqual(core, {1, 2})

    def test_periodic_body_reconstruction_ignores_bad_absolute_images(self):
        pos = np.array([[10.0, 10, 10], [39.0, 0, 0], [1.0, 0, 0]])
        mol = np.array([1, 2, 2])
        info = dict(
            np_radius_bead_centers_A=3,
            np_molecule_id_start=2,
            total_molecules=2,
            beads_per_np=2,
        )
        centers, dist = w.body_centers(pos, mol, info, np.array([40.0, 40, 40]))
        self.assertAlmostEqual(centers[0, 0], 40)
        self.assertAlmostEqual(dist[0][0, 1], 2)
        gap, which, _ = w.nearest(
            np.array([[5.0, 0, 0]]), centers, np.array([40.0, 40, 40]), 3
        )
        self.assertAlmostEqual(gap[0], 2)
        self.assertEqual(which[0], 0)

    def test_periodic_winding_includes_node_boundary_shifts(self):
        mapping = self.root / "map"
        mapping.mkdir()
        (mapping / "contours.tsv").write_text(
            "strand_id\tparent_molecule\tstatus\tfirst_node\tsecond_node\tatom_ids\n1\t1\tactive\t1\t2\t1,2\n2\t2\tactive\t2\t3\t3,4\n3\t3\tactive\t3\t1\t5,6\n"
        )
        (mapping / "nodes.tsv").write_text(
            "node_id\tkind\tactive_degree\tdangling_degree\tx_A\ty_A\tz_A\tmember_molecules\tattachment_atoms\n1\tchemical_junction\t2\t0\t0\t0\t0\t4\t1,6\n2\tchemical_junction\t2\t0\t4\t0\t0\t5\t2,3\n3\tchemical_junction\t2\t0\t8\t0\t0\t6\t4,5\n"
        )
        s = dict(
            pos=np.array(
                [
                    [0.0, 0, 0],
                    [4.0, 0, 0],
                    [4.0, 0, 0],
                    [8.0, 0, 0],
                    [8.0, 0, 0],
                    [0.0, 0, 0],
                ]
            ),
            mol=np.array([1, 1, 2, 2, 3, 3]),
            L=np.array([10.0, 10, 10]),
        )
        contours, _, weights, _, summary = w.network_map(
            mapping, s, {"composition": {"total_beads": 6}}
        )
        self.assertEqual(summary["corrected_periodic_cycle_winding_rank"], 1)
        self.assertEqual(summary["giant_component_two_core_candidate_edges"], 3)
        self.assertAlmostEqual(weights["active"].sum(), 3)

    def test_tensile_units_sign_and_known_modulus(self):
        c = self.configuration()
        n = 1000
        raw = np.zeros((n, 16))
        step = np.arange(1, n + 1) * 1000
        strain = step * 1e-4 / 1000
        raw[:, 0] = step
        raw[:, 1] = step * 5 - 2250
        raw[:, 2:4] = 300
        raw[:, 4] = strain
        raw[:, 14] = 12.3 * strain
        raw[:, 15] = raw[:, 14]
        raw[:, 7] = -raw[:, 14] / 0.101325
        raw[:, 10:13] = 100
        raw[:, 13] = 1e6
        path = self.root / "stress.dat"
        np.savetxt(path, raw)
        r = w.tensile_summary(path, "x", c, self.root)
        self.assertTrue(
            all(
                abs(f["Young_modulus_MPa"] - 12.3) < 1e-8
                for f in r["small_strain_fits"]
            )
        )
        raw[:, 14] *= -1
        np.savetxt(path, raw)
        with self.assertRaisesRegex(ValueError, "sign/unit"):
            w.tensile_summary(path, "x", c, self.root)

    def test_nearest_shell_volumes_conserve_box_and_do_not_double_count_overlap(self):
        c = self.configuration(
            volume_points=100000,
            volume_replicates=3,
            shell_edges_A=[0, 5, 10, 30],
            reference_gap_A=30,
            interface_gap_A=10,
        )
        s = dict(L=np.array([40.0, 40, 40]), low=np.zeros(3))
        centers = np.array([[15.0, 20, 20], [25.0, 20, 20]])
        v, _, per = w.volumes(s, centers, {"np_radius_bead_centers_A": 3}, c)
        self.assertAlmostEqual(v.sum(), 64000)
        np.testing.assert_allclose(per.sum(axis=0), v)
        naive = 2 * 4 * np.pi / 3 * (8**3 - 3**3)
        self.assertLess(v[1], 0.95 * naive)

    def write_dump(self, path, s, translations):
        with path.open("w") as f:
            for step, shift in enumerate(translations):
                f.write(
                    f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n4\nITEM: BOX BOUNDS pp pp pp\n0 40\n0 40\n0 40\nITEM: ATOMS id mol type x y z ix iy iz\n"
                )
                for i, xyz in enumerate((s["pos"] + shift) % 40):
                    f.write(
                        f"{i+1} {s['mol'][i]} {s['typ'][i]} {xyz[0]} {xyz[1]} {xyz[2]} 0 0 0\n"
                    )

    def test_rigid_translation_and_pbc_give_zero_relative_motion(self):
        s = dict(
            pos=np.array(
                [[39.0, 10, 10], [10.0, 10, 10], [11.0, 10, 10], [10.0, 11, 10]]
            ),
            mol=np.array([1, 2, 2, 2]),
            typ=np.array([1, 6, 6, 6]),
            L=np.array([40.0, 40, 40]),
            mass=np.ones(4),
        )
        info = dict(
            np_radius_bead_centers_A=3,
            np_molecule_id_start=2,
            total_molecules=2,
            beads_per_np=3,
            np_count=1,
            np_atom_type=6,
        )
        path = self.root / "dump"
        self.write_dump(
            path,
            s,
            [np.array([0.0, 0, 0]), np.array([2.0, 0, 0]), np.array([4.0, 0, 0])],
        )
        masks = {
            k: np.zeros(4, dtype=bool)
            for k in (
                "native_active",
                "two_core_candidate",
                "defective_strand",
                "strand_midpoint",
                "junction_attachment",
            )
        }
        c = self.configuration(dw_lags_ps=[1.0, 2.0], timestep_fs=1000.0)
        out = self.root / "dyn"
        result = w.trajectory_analysis(
            dict(path=str(path), kind="dw"), s, info, {}, masks, c, out
        )
        self.assertGreater(result["image_flag_increment_mismatch_observations"], 0)
        for row in w.read_tsv(out / "radial_dynamics.tsv"):
            if row["group"] == "matrix" and int(row["bead_origin_pairs"]):
                self.assertAlmostEqual(
                    float(row["MSD_global_COM_corrected_A2"]), 0, places=10
                )
                self.assertAlmostEqual(
                    float(row["MSD_relative_to_origin_NP_A2"]), 0, places=10
                )
        c["origin_stride"] = 2
        out = self.root / "strided"
        w.trajectory_analysis(
            dict(path=str(path), kind="dw"), s, info, {}, masks, c, out
        )
        for row in w.read_tsv(out / "NP_COM_MSD.tsv"):
            self.assertEqual(int(row["origins"]), 1)
        for row in w.read_tsv(out / "radial_dynamics.tsv"):
            self.assertEqual(int(row["origins"]), 1)
        path.write_text(path.read_text().rsplit("\n", 2)[0] + "\n")
        with self.assertRaises(ValueError):
            list(w.dump_frames(path, s))

    def test_real_metadata_mapping_and_part2_entrypoint(self):
        conf = self.root / "model.conf"
        conf.write_text(
            "system=composite\nbase=v22\noutput_dir=case\nn1=8\nm1=16\nn2=8\nfunctionality=4\nm4=0\nnp_count=1\nnp_radius=10\ndensity=.1\ntarget_density=.8\n"
        )
        subprocess.run(
            [str(ARGS.generator), "--config", str(conf)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        request = self.root / "interface.conf"
        request.write_text(
            json.dumps(
                dict(
                    data="case/data.case",
                    info="case/case.info",
                    output_dir="analysis",
                    backend_dir=str(ARGS.backend_dir),
                    plot=False,
                    volume_points=10000,
                    volume_replicates=2,
                )
            )
        )
        p = subprocess.run(
            [
                sys.executable,
                str(ROOT / "src/siliconelab_part2.py"),
                "--config",
                str(request),
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        result = json.loads((self.root / "analysis/analysis.json").read_text())
        self.assertEqual(result["status"], "complete")
        self.assertEqual(result["network"]["native_status_counts"], {"isolated": 16})
        self.assertEqual(result["network"]["corrected_periodic_cycle_winding_rank"], 0)
        rows = w.read_tsv(self.root / "analysis/equilibrium/radial_profiles.tsv")
        self.assertEqual(sum(int(r["strands_possible_sites"]) for r in rows), 32)
        self.assertTrue(
            all(float(r["native_active_strand_equivalents"]) == 0 for r in rows)
        )
        # Correlation cannot silently accept another source snapshot/metadata.
        c = w.configure(request)
        c["part2_report"] = str(self.root / "analysis/analysis.json")
        c["data"] = str(conf)
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            w.part3(c)
        c = w.configure(request)
        c["part2_report"] = str(self.root / "analysis/analysis.json")
        c["output_dir"] = str(self.root / "part3")
        self.assertEqual(w.part3(c)["status"], "complete")

    def test_affine_tensile_series_preserves_rigid_bodies_and_has_zero_residual(self):
        source = dict(
            pos=np.array(
                [
                    [20.0, 10, 10],
                    [25.0, 15, 10],
                    [10.0, 10, 10],
                    [11.0, 10, 10],
                    [10.0, 11, 10],
                ]
            ),
            mol=np.array([1, 1, 2, 2, 2]),
            typ=np.array([1, 1, 6, 6, 6]),
            L=np.array([40.0, 40, 40]),
            low=np.zeros(3),
            mass=np.ones(5),
        )
        info = dict(
            np_radius_bead_centers_A=3,
            np_molecule_id_start=2,
            total_molecules=2,
            beads_per_np=3,
            np_count=1,
            np_atom_type=6,
        )
        path = self.root / "tensile.dump"
        with path.open("w") as f:
            for step, stretch in [(0, 1.0), (1000, 1.1)]:
                pos = source["pos"].copy()
                pos[:2, 0] *= stretch
                center = pos[2:].mean(axis=0)
                pos[2:, 0] += (stretch - 1) * center[0]
                f.write(
                    f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n5\nITEM: BOX BOUNDS pp pp pp\n0 {40*stretch}\n0 40\n0 40\nITEM: ATOMS id mol type x y z ix iy iz\n"
                )
                for i, xyz in enumerate(pos):
                    f.write(
                        f"{i+1} {source['mol'][i]} {source['typ'][i]} {xyz[0]} {xyz[1]} {xyz[2]} 0 0 0\n"
                    )
        c = self.configuration(volume_points=10000, volume_replicates=2)
        weights = {
            k: (
                np.array([0.5, 0.5, 0, 0, 0])
                if k in ("all", "active", "core")
                else np.zeros(5)
            )
            for k in ("all", "active", "core", "dangling")
        }
        out = self.root / "evolution"
        r = w.tensile_interface_series(path, "x", source, info, weights, c, out)
        self.assertEqual(r["frames_checked"], 2)
        self.assertEqual(len(r["sampled_frames"]), 2)
        self.assertAlmostEqual(r["sampled_frames"][-1]["engineering_strain"], 0.1)
        for row in w.read_tsv(out / "deformation_profiles.tsv"):
            if int(row["matrix_beads"]):
                self.assertAlmostEqual(
                    float(row["minimum_image_affine_residual_MSD_A2"]), 0, places=10
                )

    def test_invalid_tensile_geometry_is_scoped_and_preserves_evidence(self):
        from unittest.mock import patch

        source = dict(
            pos=np.array([[10.0, 10, 10], [11.0, 10, 10], [10.0, 11, 10]]),
            mol=np.array([2, 2, 2]),
            typ=np.array([6, 6, 6]),
            L=np.array([40.0, 40, 40]),
            low=np.zeros(3),
            mass=np.ones(3),
        )
        info = dict(
            np_radius_bead_centers_A=3,
            np_molecule_id_start=2,
            total_molecules=2,
            beads_per_np=3,
            np_count=1,
            np_atom_type=6,
        )
        path = self.root / "bad.dump"
        path.write_text("preserved malformed rigid-body coordinates")
        pos = source["pos"].copy()
        pos[1, 0] += 0.03
        frame = (1000, pos, np.zeros((3, 3)), source["L"], source["low"])
        c = self.configuration()
        out = self.root / "invalid"
        with (
            patch.object(w, "dump_frames", return_value=iter([frame])),
            patch.object(w, "volumes") as volume,
        ):
            result = w.tensile_interface_series(path, "x", source, info, {}, c, out)
            volume.assert_not_called()
        self.assertEqual(result["status"], "invalid_rigid_geometry")
        self.assertAlmostEqual(result["maximum_NP_distance_change_A"], 0.03)
        self.assertEqual(result["invalid_frames"][0]["step"], 1000)
        self.assertEqual(result["source_sha256"], w.sha(path))
        self.assertTrue((out / "sampling.json").is_file())
        self.assertTrue((out / "rigidity_checks.tsv").is_file())
        self.assertFalse((out / "deformation_profiles.tsv").exists())
        # A user-selected tolerance is honored and saved by the normal config path.
        c = self.configuration(rigidity_tolerance_A=0.04)
        self.assertEqual(c["rigidity_tolerance_A"], 0.04)

    def test_final_np_geometry_matches_ids_after_row_permutation_and_pbc_rotation(self):
        info = dict(total_beads=5, np_atom_type=6, np_count=1, beads_per_np=4)
        source = np.array(
            [[10.0, 10, 10], [11.0, 10, 10], [10.0, 12, 10], [10.0, 10, 13]]
        )
        rotation = np.array([[1.0, 0, 0], [0.0, 0, -1], [0.0, 1, 0]])
        final = (source - 10) @ rotation.T + [39.5, 10, 10]

        def write(path, positions, order):
            rows = {1: "1 1 1 0 5 5 5 0 0 0"}
            for ident, xyz in enumerate(positions, 2):
                wrapped = xyz % 40
                image = np.floor(xyz / 40).astype(int)
                rows[ident] = f"{ident} 2 6 0 " + " ".join(map(str, [*wrapped, *image]))
            path.write_text(
                "LAMMPS data\n\n5 atoms\n2 atom types\n0 40 xlo xhi\n0 40 ylo yhi\n0 40 zlo zhi\n\nAtoms # full\n\n"
                + "\n".join(rows[i] for i in order)
                + "\n"
            )

        a, b = self.root / "source.data", self.root / "final.data"
        write(a, source, [1, 2, 3, 4, 5])
        write(b, final, [5, 1, 3, 2, 4])
        reference = np.array(d.np_geometry(a, info)[2])
        np.testing.assert_allclose(d.np_geometry(b, info)[2], reference, atol=1e-12)
        final[2, 1] += 0.2
        write(b, final, [3, 4, 1, 5, 2])
        self.assertGreater(
            float(np.max(np.abs(np.array(d.np_geometry(b, info)[2]) - reference))), 0.02
        )

    def test_added_dynamics_preserves_cold_header_and_validates_sampling(self):
        request = self.root / "dynamics.json"
        request.write_text(
            json.dumps(
                dict(
                    source_data="/case/data.final",
                    source_info="/case/case.info",
                    forcefield_input="/case/in.tensile",
                    output_dir=str(self.root / "dynamics"),
                )
            )
        )
        c = d.config(request)
        header = (
            "units real\npair_style hybrid lj/gromacs 12 15 lj/cut 7.2\nread_data old\n"
            + "neigh_modify exclude molecule/intra nanoparticles\ntimestep 5\nfix load all deform 1 x erate 1e-8\n"
        )
        text = d.input_text(c, header)
        self.assertIn("run 20000", text)
        self.assertIn("run 5000000", text)
        self.assertNotIn("fix load", text)
        self.assertNotIn("velocity ", text)
        self.assertEqual(c["dw_stride"], 20)
        self.assertEqual(c["msd_stride"], 5000)
        raw = json.loads(request.read_text())
        raw["dw_sample_ps"] = 0.003
        request.write_text(json.dumps(raw))
        with self.assertRaises(ValueError):
            d.config(request)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]] + REST)
