#!/usr/bin/env python3
"""Compute-node, disposable input preflights and rigid-body smoke checks."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess

import bulk_workflow as w


def short_input(text, steps):
    lines = []
    for line in text.splitlines():
        fields = line.split()
        if fields[:1] == ["run"]:
            if steps == 0:
                line = "run 0"
            elif int(fields[1]) > 1:
                line = "run " + str(steps)
        lines.append(line)
    return "\n".join(lines) + "\n"


def execute(folder, text, lammps):
    folder.mkdir()
    (folder / "in.validation").write_text(text)
    with (folder / "stdout.log").open("w") as out, (folder / "stderr.log").open("w") as err:
        subprocess.run(["srun", "--ntasks=1", str(lammps), "-in", "in.validation", "-log", "lammps.log"],
                       cwd=folder, stdout=out, stderr=err, check=True)
    log = (folder / "lammps.log").read_text()
    if w.has_lammps_error(log):
        raise RuntimeError("LAMMPS validation error: " + str(folder))
    return {"directory": str(folder), "warnings": sorted(set(line.strip() for line in log.splitlines() if "WARNING:" in line)),
            "executed_steps": sum(map(int, re.findall(r"Loop time of .*? for (\d+) steps", log)))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--batch-dir", type=Path, required=True)
    p.add_argument("--generator", type=Path, required=True)
    p.add_argument("--lammps", default="lmp")
    a = p.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    batch = a.batch_dir.resolve()
    plan = w.load(batch / "workflow.json")
    state = w.load(batch / "workflow-state.json")
    assert state["status"] == "prepared_inputs"
    evidence = {"status": "running", "production_inputs_modified": False, "full_case_preflights": [], "smoke_tests": []}
    w.save(batch / "preparation-result.json", evidence)
    try:
        for case in plan["cases"]:
            folder = batch / case["id"]
            validation = folder / "validation"
            validation.mkdir()
            formation = folder / "formation"
            spec = w.load(formation / "stage.json")
            reference = w.rigid_reference(spec)
            original = (formation / spec["input"]).read_text()
            text = short_input(original, 0).replace("read_data       " + spec["initial_data"], "read_data       " + str(formation / spec["initial_data"]))
            # Composite transformation normalizes this command's spacing.
            text = re.sub(r"^read_data\s+" + re.escape(spec["initial_data"]), "read_data " + str(formation / spec["initial_data"]), text, flags=re.M)
            check = validation / "formation-preflight"
            result = execute(check, text, a.lammps)
            assert result["executed_steps"] == 0
            atoms, _ = w.snapshot(check / spec["final_data"], case["expected_atoms"])
            assert w.rigidity_error({i: r[2] for i, r in atoms.items() if r[1] == 6}, reference) < 1e-6
            for axis in "xyz":
                t = (folder / ("tensile-" + axis) / "in.tensile").read_text()
                t = t.replace("../formation/" + spec["final_data"], str(check / spec["final_data"]))
                r = execute(validation / ("tensile-" + axis + "-preflight"), short_input(t, 0), a.lammps)
                assert r["executed_steps"] == 0
                result.setdefault("tensile_preflights", {})[axis] = r
            evidence["full_case_preflights"].append({"case": case["id"], **result})
            print("SUMMARY: Full-case formation and x/y/z preflights passed: " + case["id"], flush=True)
            w.save(batch / "preparation-result.json", evidence)

        scratch = batch / "validation-fixtures"
        scratch.mkdir()
        for name, lam in (("pure", None), ("native-1", 1), ("native-05", .5)):
            config = scratch / (name + ".conf")
            config.write_text(("system=coating\nmodel=v22\n" if lam is None else "system=composite\nbase=v22\n") +
                              f"output_dir={name}\nn1=24\nm1=24\nn2=16\nfunctionality=4\nm4=0\ndensity=.1\ntarget_density=.8\n" +
                              ("" if lam is None else f"np_count=2\nnp_radius=10\nnp_interaction=native_pdms\nnp_lambda={lam}\n"))
            subprocess.run([str(a.generator), "--config", str(config)], check=True)
            package = scratch / name
            original = (package / ("in." + name)).read_text()
            output = "data." + name + ".validation_final"
            text = short_input(original, 20) + "\nwrite_data " + output + " nocoeff\n"
            text = re.sub(r"^read_data\s+data\." + re.escape(name), "read_data " + str(package / ("data." + name)), text, flags=re.M)
            run = execute(scratch / (name + "-formation-smoke"), text, a.lammps)
            count = int(next(line.split()[0] for line in (package / ("data." + name)).read_text().splitlines() if line.endswith(" atoms")))
            spec = {"lambda": lam, "package": str(package), "initial_data": "data." + name, "expected_atoms": count}
            # The production identity contract is intentionally 32x429. Small fixtures use their own body sizes.
            reference = {}
            if lam is not None:
                initial, _ = w.snapshot(package / ("data." + name), count)
                groups = {}
                for i, (mol, typ, xyz) in initial.items():
                    if typ == 6:
                        groups.setdefault(mol, []).append((i, xyz))
                import math
                reference = {mol: {"ids": [r[0] for r in rows], "anchors": [r[0] for r in rows[:4]],
                                   "distances": [[math.dist(a[1], r[1]) for r in rows] for a in rows[:4]]}
                             for mol, rows in groups.items()}
            atoms, _ = w.snapshot(Path(run["directory"]) / output, count)
            assert w.rigidity_error({i: r[2] for i, r in atoms.items() if r[1] == 6}, reference) < 1e-6
            for axis in "xyz":
                t = w.tensile_input(original, str(Path(run["directory"]) / output), axis, lam is not None, 2e-8, 1)
                d = scratch / (name + "-" + axis + "-smoke")
                r = execute(d, short_input(t, 20), a.lammps)
                atoms, _ = w.snapshot(d / "data.tensile_final", count)
                error = w.rigidity_error({i: row[2] for i, row in atoms.items() if row[1] == 6}, reference)
                assert error < 1e-6
                run.setdefault("tensile_smokes", {})[axis] = {**r, "maximum_NP_distance_error_A": error}
            evidence["smoke_tests"].append({"case": name, **run})
            w.save(batch / "preparation-result.json", evidence)
        state.update(status="validated", validated_at=w.now())
        w.save(batch / "workflow-state.json", state)
        evidence.update(status="complete", validated_at=w.now())
    except Exception as e:
        evidence.update(status="failed", error=str(e))
        raise
    finally:
        w.save(batch / "preparation-result.json", evidence)
    print("SUMMARY: All seven production packages and disposable formation/tensile smoke tests verified; production inputs unchanged.", flush=True)


if __name__ == "__main__":
    main()
