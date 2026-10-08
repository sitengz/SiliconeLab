#!/usr/bin/env python3
"""Compute-node preparation and guarded, resumable Slurm film DAG."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "workflows/films.json"
STATE = ROOT / "runs/film-workflow.json"
MODULE = "lammps/20230802.2-py310-openmpi4-ezoqd7f"


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def load(path):
    return json.loads(path.read_text())


def now():
    return datetime.now(timezone.utc).isoformat()


def rootpath(p):
    return ROOT / p


def data_header(path):
    counts = {}
    bounds = {}
    with path.open() as f:
        for line in f:
            words = line.split()
            if words and words[0] == "Atoms":
                break
            if len(words) == 2 and words[1] in (
                "atoms",
                "bonds",
                "angles",
                "dihedrals",
            ):
                counts[words[1]] = int(words[0])
            if len(words) == 4 and words[2] in ("xlo", "ylo", "zlo"):
                lo, hi = map(float, words[:2])
                assert hi > lo
                bounds[words[2][0]] = [lo, hi]
    if set(bounds) != set("xyz") or counts.get("atoms", 0) < 1:
        raise RuntimeError("Invalid or missing data header: " + str(path))
    return counts, bounds


def cold_header(native, source):
    """Select the final 300 K coefficients, preserving native force-field styles."""
    lines = native.read_text().splitlines()
    last = {}
    ordered = []
    styles = (
        "units",
        "atom_style",
        "bond_style",
        "angle_style",
        "dihedral_style",
        "special_bonds",
        "pair_style",
        "comm_modify",
        "neighbor",
        "neigh_modify",
    )
    for line in lines:
        w = line.split()
        if not w:
            continue
        if w[0] in styles:
            key = (w[0],)
        elif w[0] in (
            "bond_coeff",
            "angle_coeff",
            "dihedral_coeff",
            "pair_coeff",
            "mass",
        ):
            key = tuple(w[:3] if w[0] == "pair_coeff" else w[:2])
        else:
            continue
        if key not in last:
            ordered.append(key)
        last[key] = line
    if not any(k[0] == "pair_coeff" for k in ordered):
        raise RuntimeError("No native pair matrix")
    before = [last[k] for k in ordered if k[0] in styles]
    after = [last[k] for k in ordered if k[0] not in styles]
    # read_data must follow styles and precede coefficients. No bulk network is cut.
    return (
        "\n".join(
            before[:1]
            + ["boundary p p f"]
            + before[1:]
            + ["read_data " + source]
            + after
            + [
                "timestep 5",
                "thermo 1000",
                "thermo_style custom step time temp pe pxx pyy pzz lx ly lz",
                "thermo_modify lost error flush yes",
                "restart 1000000 restart.stage.1 restart.stage.2",
            ]
        )
        + "\n"
    )


def guards(expand):
    t = f"change_box all z delta -{expand} {expand} units box\n" if expand else ""
    for side in ("zlo", "zhi"):
        t += f"fix {side}_wall all wall/lj126 {side} EDGE 1.012878448 6.445843658 7.235214876 units box\nfix_modify {side}_wall virial no\n"
    return t


def pressure_output(filename, label):
    return (
        "\n".join(
            [
                "variable t equal time",
                "variable T equal temp",
                "variable e equal pe",
                "variable x equal pxx",
                "variable y equal pyy",
                "variable z equal pzz",
                "variable a equal lx",
                "variable b equal ly",
                "variable c equal lz",
                "variable wl equal f_zlo_wall[1]",
                "variable wh equal f_zhi_wall[1]",
                f'fix {label} all print 1000 "${{t}} ${{T}} ${{e}} ${{x}} ${{y}} ${{z}} ${{a}} ${{b}} ${{c}} ${{wl}} ${{wh}}" file {filename} screen no title "# time_fs temp_K pe_kcal_per_mol pxx_atm pyy_atm pzz_atm lx_A ly_A lz_A wall_lo_force wall_hi_force"',
            ]
        )
        + "\n"
    )


def build_surface(path, header, expand):
    text = header + guards(expand)
    for phase in ("equil", "production"):
        text += "reset_timestep 0 time 0.0\nfix integrate all nvt temp 300 300 50\n"
        text += pressure_output("pressure." + phase + ".dat", "pressure")
        text += f"dump geom all custom 100000 dump.{phase}.lammpstrj id mol type x y z ix iy iz\ndump_modify geom first yes sort id\nrun 10000000\nunfix pressure\nunfix integrate\nundump geom\nwrite_data data.surface_{phase} nocoeff\n"
    text += 'print "SILICONELAB_STAGE_COMPLETE surface"\n'
    path.write_text(text)


def build_dynamics(path, header):
    text = header + guards(0)
    text += "fix integrate all nvt temp 300 300 50\nrun 1000000\nwrite_data data.dynamics_equil nocoeff\nreset_timestep 0 time 0.0\n"
    text += "dump dw all custom 20 dump.debye_waller.lammpstrj id mol type x y z ix iy iz\ndump_modify dw first yes sort id\nrun 20000\nundump dw\nreset_timestep 0 time 0.0\n"
    text += 'dump msd all custom 1000 dump.layer_dynamics.lammpstrj id mol type x y z ix iy iz\ndump_modify msd first yes sort id\nrun 1000000\ndump_modify msd every 5000 first no\nrun 4000000\nundump msd\nunfix integrate\nwrite_data data.dynamics_final nocoeff\nprint "SILICONELAB_STAGE_COMPLETE dynamics"\n'
    path.write_text(text)


def build_tensile(path, header, axis, h):
    other = "y" if axis == "x" else "x"
    text = header + guards(0)
    # Fixed free-surface cell z; lateral transverse pressure control only.
    text += f"fix eq all npt temp 300 300 50 x 1 1 500 y 1 1 500 couple xy\nrun 1000000\nunfix eq\nrun 0\nvariable Lx0 equal $(lx:%.12g)\nvariable Ly0 equal $(ly:%.12g)\nvariable Lz0 equal $(lz:%.12g)\n"
    text += 'print "Lx0_A Ly0_A Lz0_A" file reference_box.dat screen no\nprint "${Lx0} ${Ly0} ${Lz0}" append reference_box.dat screen no\nreset_timestep 0 time 0.0\n'
    text += f"variable t equal time\nvariable T equal temp\nvariable ex equal lx/v_Lx0-1\nvariable ey equal ly/v_Ly0-1\nvariable ez equal lz/v_Lz0-1\nvariable x equal pxx\nvariable y equal pyy\nvariable z equal pzz\nvariable a equal lx\nvariable b equal ly\nvariable c equal lz\nvariable v equal vol\n"
    text += f"fix integrate all npt temp 300 300 50 {other} 1 1 500\nfix deform all deform 1 {axis} erate 2e-9 units box remap x\n"
    text += "fix record all ave/time 100 10 1000 v_t v_T v_ex v_ey v_ez v_x v_y v_z v_a v_b v_c v_v file stress_strain.dat\n"
    text += 'dump geom all custom 1000000 dump.tensile.lammpstrj id mol type x y z ix iy iz\ndump_modify geom first yes sort id\nrun 20000000\nunfix record\nunfix deform\nunfix integrate\nundump geom\nwrite_data data.tensile_final nocoeff\nprint "SILICONELAB_STAGE_COMPLETE tensile"\n'
    path.write_text(text)


def add_stage(state, case, stage, directory, input_name, final, dependencies, required):
    key = case + ":" + stage
    item = {
        "case": case,
        "stage": stage,
        "directory": str(directory.relative_to(ROOT)),
        "input": input_name,
        "final_data": final,
        "dependencies": dependencies,
        "required_outputs": required,
        "state": "prepared",
        "job_id": None,
        "expected_atoms": state["cases"][case]["atoms"],
    }
    item["input_sha256"] = digest(directory / input_name)
    save(
        directory / "stage-metadata.json",
        dict(
            item,
            geometry="film",
            timestep_fs=5,
            source_info=state["cases"][case]["source_info"],
            nominal_material_thickness_A=state["cases"][case]["nominal_thickness_A"],
            surface_definition="component mass-density half-maximum crossing; never simulation box faces",
            sampling=load(PLAN)["sampling"],
        ),
    )
    state["stages"][key] = item
    script = directory / "submit.sbatch"
    script.write_text(
        "#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks=96\n#SBATCH --mem=200G\n#SBATCH --time=48:00:00\n#SBATCH --account=wxia\n#SBATCH --partition=nova\nset -eo pipefail\nsource /etc/profile\nsource /work/wxia/siteng/log/tools/nova-logging.sh\nmodule purge\nmodule load "
        + MODULE
        + '\nexport OMP_NUM_THREADS=1\nnova_run "Execute and verify film workflow stage" --summary "Compute-node stage records inputs, parent checks, LAMMPS completion and output verification in stage-result.json." python3 '
        + shlex.quote(str(ROOT / "scripts/film_workflow.py"))
        + " run "
        + shlex.quote(key)
        + "\n"
    )


def prepare():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Preparation must run on a compute node")
    if STATE.exists():
        raise RuntimeError("Preparation already recorded; use saved state")
    plan = load(PLAN)
    state = {
        "schema_version": 1,
        "created_at": now(),
        "revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "cases": {},
        "stages": {},
    }
    for c in plan["cases"]:
        case = c["id"]
        directory = rootpath(c["run_directory"])
        bulk = rootpath(c["bulk_run"])
        bc = bulk.name
        parent = load(bulk / "submission.json")
        # sacct is authoritative; old metadata is not itself a completion gate.
        accounting = subprocess.check_output(
            [
                "sacct",
                "-j",
                str(c["bulk_job_id"]),
                "-n",
                "-P",
                "--format=JobIDRaw,State,ExitCode",
            ],
            text=True,
        )
        if f"{c['bulk_job_id']}|COMPLETED|0:0" not in accounting:
            raise RuntimeError("Bulk prerequisite failed: " + bc)
        bulk_data = bulk / f"data.{bc}.npt_eq"
        counts, bounds = data_header(bulk_data)
        if counts["atoms"] != c["atoms"]:
            raise RuntimeError("Bulk atom count mismatch")
        subprocess.run(
            [
                str(ROOT / "build/siliconelab_generator"),
                "--config",
                str(rootpath(c["config"])),
            ],
            check=True,
        )
        native = directory / ("in." + case + (".film" if c["system"] == "oil" else ""))
        if c["system"] == "oil":
            shutil.copy2(bulk_data, directory / f"data.{case}.npt_eq")
            shutil.copy2(bulk / (bc + ".info"), directory / "bulk-source.info")
            # Keep native conversion and full 50+50 ns. Add geometry dumps and marker.
            t = native.read_text()
            for phase in ("equil", "production"):
                token = (
                    "fix             "
                    + ("equil_output" if phase == "equil" else "energy_output")
                    + " all print"
                )
                dump = f"dump geometry all custom 100000 dump.surface_{phase}.lammpstrj id mol type x y z ix iy iz\ndump_modify geometry first yes sort id\n"
                t = t.replace(token, dump + token, 1)
                endfix = "unfix           " + (
                    "equil_output" if phase == "equil" else "energy_output"
                )
                t = t.replace(endfix, endfix + "\nundump geometry", 1)
            t += '\nprint "SILICONELAB_STAGE_COMPLETE film"\n'
            native.write_text(t)
            final = "data." + case + ".film_final"
            required = [
                final,
                "energy." + case + ".film.dat",
                "energy." + case + ".film_eq.dat",
                "dump.surface_production.lammpstrj",
            ]
        else:
            native.write_text(
                native.read_text() + '\nprint "SILICONELAB_STAGE_COMPLETE film"\n'
            )
            final = "data." + case + ".npt_eq"
            required = [final, "dump.msd.lammpstrj"]
        source_info = load(directory / (case + ".info"))
        state["cases"][case] = {
            "atoms": c["atoms"],
            "system": c["system"],
            "bulk_case": bc,
            "source_info": str((directory / (case + ".info")).relative_to(ROOT)),
            "bulk_data_sha256": digest(bulk_data),
            "nominal_thickness_A": c["thickness_A"],
            "config_sha256": digest(rootpath(c["config"])),
        }
        add_stage(state, case, "film", directory, native.name, final, [], required)
        if c["system"] == "oil":
            surface_key = case + ":film"
            surface_data = directory / final
        else:
            surf = directory / "surface"
            surf.mkdir()
            surface_key = case + ":surface"
            surface_data = surf / "data.surface_production"
            header = cold_header(native, "../" + final)
            build_surface(surf / "in.surface", header, 50)
            add_stage(
                state,
                case,
                "surface",
                surf,
                "in.surface",
                "data.surface_production",
                [case + ":film"],
                [
                    "data.surface_production",
                    "pressure.equil.dat",
                    "pressure.production.dat",
                    "dump.production.lammpstrj",
                ],
            )
        # Headers use a relative source; each independent property stage reads the same equilibrated surface state.
        for stage in ["dynamics"] + (
            [] if c["system"] == "oil" else ["tensile-x", "tensile-y"]
        ):
            sub = directory / stage
            sub.mkdir()
            rel = os.path.relpath(surface_data, sub)
            header = cold_header(native, rel)
            if stage == "dynamics":
                build_dynamics(sub / "in.dynamics", header)
                fname = "data.dynamics_final"
                req = [
                    fname,
                    "data.dynamics_equil",
                    "dump.debye_waller.lammpstrj",
                    "dump.layer_dynamics.lammpstrj",
                ]
                input_name = "in.dynamics"
            else:
                build_tensile(sub / "in.tensile", header, stage[-1], c["thickness_A"])
                fname = "data.tensile_final"
                req = [
                    fname,
                    "stress_strain.dat",
                    "reference_box.dat",
                    "dump.tensile.lammpstrj",
                ]
                input_name = "in.tensile"
            add_stage(state, case, stage, sub, input_name, fname, [surface_key], req)
        check = directory / "preflight"
        check.mkdir()
        validation = (
            cold_header(native, str(directory / ("data." + case)))
            + guards(50)
            + "run 0\n"
        )
        (check / "in.forcefield").write_text(validation)
        subprocess.run(
            [
                "srun",
                "--ntasks=1",
                "lmp",
                "-in",
                "in.forcefield",
                "-log",
                "lammps.log",
                "-screen",
                "screen.out",
            ],
            cwd=check,
            check=True,
        )
        if "ERROR" in (check / "lammps.log").read_text():
            raise RuntimeError("Force-field preflight failed")
        state["cases"][case]["forcefield_preflight"] = "passed"
        save(STATE, state)
    # Prevent use of a partially prepared graph after a preparation failure.
    state["preparation_complete"] = True
    save(STATE, state)
    print(
        "SUMMARY: Four full-size film packages and 17 stages prepared; input hashes and bulk prerequisites recorded."
    )


def check_stage(item):
    directory = rootpath(item["directory"])
    result = directory / "stage-result.json"
    if not result.exists() or load(result).get("status") != "complete":
        raise RuntimeError("Missing completion contract: " + str(result))
    for name in item["required_outputs"]:
        p = directory / name
        if not p.is_file() or p.stat().st_size < 1:
            raise RuntimeError("Missing output: " + str(p))
    counts, _ = data_header(directory / item["final_data"])
    if counts["atoms"] != item["expected_atoms"]:
        raise RuntimeError("Final atom count mismatch")


def run_stage(key):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("MD must run on a compute node")
    state = load(STATE)
    item = state["stages"][key]
    directory = rootpath(item["directory"])
    for parent in item["dependencies"]:
        pitem = state["stages"][parent]
        check_stage(pitem)
        if (
            digest(rootpath(pitem["directory"]) / pitem["final_data"])
            != load(rootpath(pitem["directory"]) / "stage-result.json")[
                "final_data_sha256"
            ]
        ):
            raise RuntimeError("Parent snapshot changed after completion")
    if digest(directory / item["input"]) != item["input_sha256"]:
        raise RuntimeError("Scientific input changed after preparation")
    if (directory / "stage-result.json").exists():
        raise RuntimeError("Refusing to overwrite an existing stage result")
    result = {
        "status": "running",
        "key": key,
        "job_id": os.environ["SLURM_JOB_ID"],
        "started_at": now(),
        "input_sha256": item["input_sha256"],
    }
    path = directory / "stage-result.json"
    try:
        # srun MPI dependencies come exclusively from the validated LAMMPS module.
        subprocess.run(
            ["srun", "--ntasks=96", "lmp", "-in", item["input"], "-log", "lammps.log"],
            cwd=directory,
            check=True,
        )
        log = (directory / "lammps.log").read_text()
        if (
            "ERROR" in log
            or "SILICONELAB_STAGE_COMPLETE " + item["stage"].split("-")[0] not in log
        ):
            raise RuntimeError("LAMMPS completion marker missing or error present")
        for name in item["required_outputs"]:
            p = directory / name
            if not p.is_file() or not p.stat().st_size:
                raise RuntimeError("Required output missing: " + str(p))
        counts, bounds = data_header(directory / item["final_data"])
        if counts["atoms"] != item["expected_atoms"]:
            raise RuntimeError("Final atom count mismatch")
        result.update(
            status="complete",
            final_counts=counts,
            final_bounds=bounds,
            final_data_sha256=digest(directory / item["final_data"]),
            finished_at=now(),
        )
        validate_trajectories(item, directory, result)
        save(path, result)
    except Exception as e:
        result.update(status="failed", error=str(e), finished_at=now())
        save(path, result)
        raise


def validate_trajectories(item, directory, result):
    checks = {}
    for name in item["required_outputs"]:
        if not name.endswith(".lammpstrj"):
            continue
        p = directory / name
        with p.open("rb") as f:
            first = [f.readline().decode().strip() for _ in range(9)]
            f.seek(max(0, p.stat().st_size - 32 * 1024 * 1024))
            tail = f.read().decode(errors="strict")
        steps = re.findall(r"ITEM: TIMESTEP\n([0-9]+)\n", tail)
        if (
            not steps
            or first[2] != "ITEM: NUMBER OF ATOMS"
            or int(first[3]) != item["expected_atoms"]
        ):
            raise RuntimeError("Incomplete trajectory " + name)
        if first[4] != "ITEM: BOX BOUNDS pp pp ff":
            raise RuntimeError("Wrong film trajectory boundaries: " + first[4])
        if first[8] != "ITEM: ATOMS id mol type x y z ix iy iz":
            raise RuntimeError("Unexpected trajectory columns")
        final = int(steps[-1])
        expected = (
            (20000 if "debye" in name else 5000000)
            if item["stage"] == "dynamics"
            else (
                20000000
                if item["stage"].startswith("tensile")
                else (10000000 if "production" in name else None)
            )
        )
        if expected and final != expected:
            raise RuntimeError(
                "Trajectory duration incomplete: " + name + " " + str(final)
            )
        checks[name] = {
            "initial_timestep": int(first[1]),
            "last_timestep": final,
            "size_bytes": p.stat().st_size,
            "columns": first[8],
        }
    result["trajectory_checks"] = checks


def submit():
    # Atomic accepted-job persistence and Slurm history checks prevent duplicates after interruption.
    state = load(STATE)
    if not state.get("preparation_complete"):
        raise RuntimeError("Preparation incomplete")
    for key, item in state["stages"].items():
        if item.get("job_id"):
            continue
        queue = subprocess.check_output(
            ["squeue", "--array", "-u", os.environ["USER"], "-h", "-o", "%i|%E|%j"],
            text=True,
        ).splitlines()
        independent = [r for r in queue if r.split("|")[1] in ("", "(null)", "N/A")]
        dependencies = [str(state["stages"][p]["job_id"]) for p in item["dependencies"]]
        if not dependencies and len(independent) >= 8:
            raise RuntimeError(
                "Eight independent jobs already queued/running; resume later"
            )
        directory = rootpath(item["directory"])
        name = ("SLfilm-" + key).replace(":", "-")
        if any(r.split("|")[-1] == name for r in queue):
            raise RuntimeError(
                "Unrecorded matching job found; reconcile before resubmission: " + name
            )
        historical = subprocess.check_output(
            [
                "sacct",
                "-S",
                state["created_at"][:10],
                "-n",
                "-P",
                "--name",
                name,
                "--format=JobIDRaw,JobName,State",
            ],
            text=True,
        ).strip()
        if historical:
            raise RuntimeError(
                "Historical matching job found; reconcile before resubmission: " + name
            )
        cmd = [
            "sbatch",
            "--parsable",
            "--job-name=" + name,
            "--chdir=" + str(directory),
            "--output=" + str(directory / "%j.out"),
            "--error=" + str(directory / "%j.err"),
        ]
        if dependencies:
            cmd += [
                "--dependency=afterok:" + ":".join(dependencies),
                "--kill-on-invalid-dep=yes",
            ]
        cmd += ["submit.sbatch"]
        job = (
            subprocess.check_output(cmd, cwd=directory, text=True).strip().split(";")[0]
        )
        if not job.isdigit():
            raise RuntimeError("Unexpected sbatch response")
        item.update(
            job_id=job,
            state="SUBMITTED",
            submitted_at=now(),
            slurm_stdout=str(directory / (job + ".out")),
            slurm_stderr=str(directory / (job + ".err")),
            submission_command=cmd,
        )
        save(directory / "submission.json", item)
        save(STATE, state)
        print(
            "SUMMARY: Submitted "
            + key
            + " job "
            + job
            + "; parents "
            + (",".join(dependencies) or "verified completed bulk/preparation"),
            flush=True,
        )


def monitor():
    state = load(STATE)
    for key, item in state["stages"].items():
        job = item.get("job_id")
        if not job:
            continue
        raw = subprocess.check_output(
            [
                "sacct",
                "-j",
                job,
                "-n",
                "-P",
                "--format=JobIDRaw,State,ExitCode,Elapsed",
            ],
            text=True,
        )
        row = next(
            (r.split("|") for r in raw.splitlines() if r.startswith(job + "|")), None
        )
        if not row:
            continue
        item.update(
            state=row[1], exit_code=row[2], elapsed=row[3], last_checked_at=now()
        )
        if row[1] == "COMPLETED" and row[2] == "0:0":
            try:
                check_stage(item)
                item["outputs_verified"] = True
            except Exception as e:
                item.update(outputs_verified=False, verification_error=str(e))
        print(
            "SUMMARY: "
            + key
            + " job "
            + job
            + " "
            + row[1]
            + " exit "
            + row[2]
            + " outputs_verified="
            + str(item.get("outputs_verified", False))
        )
    for key, item in state.get("analysis_jobs", {}).items():
        job = item["job_id"]
        raw = subprocess.check_output(
            [
                "sacct",
                "-j",
                job,
                "-n",
                "-P",
                "--format=JobIDRaw,State,ExitCode,Elapsed",
            ],
            text=True,
        )
        row = next(
            (r.split("|") for r in raw.splitlines() if r.startswith(job + "|")), None
        )
        if not row:
            continue
        item.update(
            state=row[1], exit_code=row[2], elapsed=row[3], last_checked_at=now()
        )
        report = rootpath(item["output_directory"]) / "analysis.json"
        if row[1] == "COMPLETED" and row[2] == "0:0":
            if report.exists():
                data = load(report)
                item["outputs_verified"] = data.get("status") == "complete"
                item["observable_statuses"] = {
                    r["analysis"]: r["status"] for r in data.get("analyses", [])
                }
            else:
                item.update(
                    outputs_verified=False, verification_error="Missing analysis.json"
                )
        print(
            "SUMMARY: "
            + key
            + " job "
            + job
            + " "
            + row[1]
            + " exit "
            + row[2]
            + " outputs_verified="
            + str(item.get("outputs_verified", False))
        )
    save(STATE, state)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["prepare", "submit", "monitor", "run"])
    p.add_argument("key", nargs="?")
    a = p.parse_args()
    {
        "prepare": prepare,
        "submit": submit,
        "monitor": monitor,
        "run": lambda: run_stage(a.key),
    }[a.action]()
