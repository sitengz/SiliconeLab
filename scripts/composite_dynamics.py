#!/usr/bin/env python3
"""Prepare and verify one additional rigid-composite NVT dynamics allocation.

All scientific and resource controls are saved in a JSON configuration. Preparation
writes text only; LAMMPS preflights, MD and complete trajectory checks run in Slurm.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

DEFAULTS = dict(
    temperature_K=300.0,
    timestep_fs=5.0,
    thermostat_damping_fs=50.0,
    dw_duration_ps=100.0,
    dw_sample_ps=0.1,
    msd_duration_ns=25.0,
    msd_sample_ps=25.0,
    checkpoint_steps=100000,
    thermo_steps=1000,
    tasks=96,
    memory="200G",
    time_limit="48:00:00",
    account="wxia",
    partition="nova",
    independent_job_limit=8,
    lammps_module="lammps/20230802.2-py310-openmpi4-ezoqd7f",
    logging_helper="/work/wxia/siteng/log/tools/nova-logging.sh",
    lammps="lmp",
    final_rigidity_tolerance_A=1e-6,
)


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def save(p, v):
    p = Path(p)
    t = p.with_name(p.name + ".tmp")
    t.write_text(json.dumps(v, indent=2, allow_nan=False) + "\n")
    t.replace(p)


def config(path):
    raw = json.loads(Path(path).read_text())
    required = {"source_data", "source_info", "forcefield_input", "output_dir"}
    if set(raw) - set(DEFAULTS) - required:
        raise ValueError("Unknown dynamics configuration key")
    if not required.issubset(raw):
        raise ValueError("Missing source/output settings")
    c = {**DEFAULTS, **raw}
    for key in required:
        c[key] = str((Path(path).resolve().parent / c[key]).resolve())
    for k in (
        "temperature_K",
        "timestep_fs",
        "thermostat_damping_fs",
        "dw_duration_ps",
        "dw_sample_ps",
        "msd_duration_ns",
        "msd_sample_ps",
        "final_rigidity_tolerance_A",
    ):
        if (
            isinstance(c[k], bool)
            or not isinstance(c[k], (float, int))
            or not math.isfinite(c[k])
            or c[k] <= 0
        ):
            raise ValueError("Invalid " + k)
    for k in ("checkpoint_steps", "thermo_steps", "tasks", "independent_job_limit"):
        if type(c[k]) is not int or c[k] <= 0:
            raise ValueError("Invalid " + k)
    if c["independent_job_limit"] > 8:
        raise ValueError("Maximum authorized concurrency is eight")
    for key, duration in [
        ("dw", c["dw_duration_ps"] * 1000),
        ("msd", c["msd_duration_ns"] * 1e6),
    ]:
        n = duration / c["timestep_fs"]
        stride = c[key + "_sample_ps"] * 1000 / c["timestep_fs"]
        if not n.is_integer() or not stride.is_integer() or int(n) % int(stride):
            raise ValueError("Duration/sample interval must fit complete timesteps")
        c[key + "_steps"] = int(n)
        c[key + "_stride"] = int(stride)
    if c["dw_steps"] % c["msd_stride"]:
        raise ValueError("DW phase end must align with MSD dump interval")
    return c


def input_text(c, source):
    lines = source.splitlines()
    header = []
    for line in lines:
        w = line.split("#", 1)[0].split()
        if w[:1] == ["timestep"]:
            break
        if w and w[0] in ("run", "fix", "velocity", "dump", "restart"):
            raise ValueError(
                "Supply an established cold input with a clean force-field header"
            )
        if w[:1] == ["read_data"]:
            line = "read_data " + str(Path(c["source_data"]).resolve())
        header.append(line)
    else:
        raise ValueError("Cold input has no timestep boundary")
    h = "\n".join(header)
    if (
        "lj/gromacs" not in h
        or "neigh_modify exclude molecule/intra nanoparticles" not in h
    ):
        raise ValueError("Expected validated cold rigid-composite header")
    T, damp = c["temperature_K"], c["thermostat_damping_fs"]
    return h + f"""
# Fixed-volume continuation of the undeformed NPT snapshot; velocities retained.
timestep {c['timestep_fs']}
reset_timestep 0
thermo {c['thermo_steps']}
variable cvx equal vcm(all,x)
variable cvy equal vcm(all,y)
variable cvz equal vcm(all,z)
thermo_style custom step time c_precursor_temp c_composite_temp density lx ly lz c_composite_press[1] c_composite_press[2] c_composite_press[3] etotal epair ebond eangle edihed v_cvx v_cvy v_cvz
thermo_modify lost error flush yes
fix np_integrate nanoparticles rigid/nvt/small molecule temp {T} {T} {damp} reinit no
fix matrix_integrate precursor nvt temp {T} {T} {damp}
run 0
restart {c['checkpoint_steps']} restart.dynamics.1 restart.dynamics.2
dump dw all custom {c['dw_stride']} dump.dw.lammpstrj id mol type x y z ix iy iz
dump_modify dw first yes sort id format float %.12g
run {c['dw_steps']}
undump dw
write_restart restart.after_dw
dump msd all custom {c['msd_stride']} dump.msd.lammpstrj id mol type x y z ix iy iz
dump_modify msd first yes sort id format float %.12g
run {c['msd_steps']}
undump msd
write_data data.dynamics_final nocoeff
write_restart restart.dynamics_final
print "SILICONELAB_COMPOSITE_DYNAMICS_COMPLETE" file completion.txt screen yes
unfix matrix_integrate
unfix np_integrate
"""


def prepare(path):
    c = config(path)
    out = Path(c["output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    source = Path(c["forcefield_input"]).read_text()
    text = input_text(c, source)
    (out / "in.dynamics").write_text(text)
    save(out / "model.conf", c)
    shutil.copy2(__file__, out / "runtime.py")
    info = json.loads(Path(c["source_info"]).read_text())
    contract = dict(
        config=c,
        input_sha256=sha(out / "in.dynamics"),
        runtime_sha256=sha(out / "runtime.py"),
        sources={
            k: {"path": c[k], "sha256": sha(c[k])}
            for k in ("source_data", "source_info", "forcefield_input")
        },
        expected_atoms=info["total_beads"],
        expected_frames={
            k: c[k + "_steps"] // c[k + "_stride"] + 1 for k in ("dw", "msd")
        },
    )
    save(out / "case.json", contract)
    script = f"""#!/bin/bash
#SBATCH --job-name=SL-NP-interface-dynamics
#SBATCH --account={c['account']}
#SBATCH --partition={c['partition']}
#SBATCH --nodes=1
#SBATCH --ntasks={c['tasks']}
#SBATCH --mem={c['memory']}
#SBATCH --time={c['time_limit']}
#SBATCH --chdir={out}
#SBATCH --output={out}/%j.out
#SBATCH --error={out}/%j.err
set -eo pipefail
source /etc/profile
source {shlex.quote(c['logging_helper'])}
nova_run 'Load established rigid-composite simulation module' --summary 'Loaded validated LAMMPS and matching MPI dependencies' module load {shlex.quote(c['lammps_module'])}
export OMP_NUM_THREADS=1
nova_run 'Run authorized NP interface dynamics and verify output contracts' --summary 'Preflight and full MD output remain in the dynamics folder; original formation and tensile files preserved' python3 {shlex.quote(str(out / 'runtime.py'))} run --config {shlex.quote(str(out / 'model.conf'))}
"""
    (out / "submit.sbatch").write_text(script)
    print(
        "SUMMARY: Prepared one new dynamics folder with saved protocol and source hashes."
    )


def frames(path, expected_atoms, expected):
    steps = []
    with Path(path).open() as f:
        while True:
            s = f.readline()
            if not s:
                break
            if s.strip() != "ITEM: TIMESTEP":
                raise ValueError("Malformed/truncated dump")
            steps.append(int(f.readline()))
            if (
                f.readline().strip() != "ITEM: NUMBER OF ATOMS"
                or int(f.readline()) != expected_atoms
            ):
                raise ValueError("Atom count changed")
            if f.readline().strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError("Bulk periodic box required")
            for _ in range(3):
                lo, hi = map(float, f.readline().split())
                assert math.isfinite(lo) and math.isfinite(hi) and hi > lo
            if f.readline().split()[2:] != [
                "id",
                "mol",
                "type",
                "x",
                "y",
                "z",
                "ix",
                "iy",
                "iz",
            ]:
                raise ValueError("Unexpected dump columns")
            for i in range(1, expected_atoms + 1):
                row = f.readline().split()
                if (
                    len(row) != 9
                    or int(row[0]) != i
                    or not all(math.isfinite(float(x)) for x in row[3:6])
                ):
                    raise ValueError("Incomplete/nonfinite atom record")
    if steps != list(expected):
        raise ValueError("Trajectory duration/sampling contract failed")
    return dict(
        frames=len(steps), first_step=steps[0], last_step=steps[-1], sha256=sha(path)
    )


def np_geometry(path, info):
    box = {}
    groups = {}
    section = None
    count = None
    atoms = 0
    with Path(path).open() as f:
        for line in f:
            w = line.split("#", 1)[0].split()
            if not w:
                continue
            if len(w) == 2 and w[1] == "atoms":
                count = int(w[0])
            if len(w) == 4 and w[2] in ("xlo", "ylo", "zlo"):
                box[w[2][0]] = float(w[1]) - float(w[0])
            if w[0] == "Atoms":
                section = "atoms"
                continue
            if not w[0].lstrip("-").isdigit():
                section = None
                continue
            if section == "atoms":
                atoms += 1
                if int(w[2]) == info["np_atom_type"]:
                    xyz = [
                        float(w[4 + j])
                        + (int(w[7 + j]) * box[a] if len(w) == 10 else 0)
                        for j, a in enumerate("xyz")
                    ]
                    groups.setdefault(int(w[1]), []).append((int(w[0]), xyz))
    if (
        count != info["total_beads"]
        or atoms != count
        or len(groups) != info["np_count"]
        or any(len(g) != info["beads_per_np"] for g in groups.values())
    ):
        raise ValueError("Final snapshot/NP identity contract failed")
    # LAMMPS write_data emits processor/spatial atom order, which can change.
    # Anchor identities and comparison positions must follow stable atom IDs.
    for rows in groups.values():
        rows.sort(key=lambda row: row[0])
    return {
        mol: [math.dist(anchor[1], r[1]) for anchor in rows[:4] for r in rows]
        for mol, rows in groups.items()
    }


def run(path):
    # model.conf is already resolved and includes derived step counts.
    c = json.loads(Path(path).read_text())
    out = Path(c["output_dir"])
    spec = json.loads((out / "case.json").read_text())
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    if (out / "stage-result.json").exists() or (out / "lammps.log").exists():
        raise ValueError("Existing evidence; no automatic MD resubmission")
    result = {"status": "running", "job_id": os.environ["SLURM_JOB_ID"], "protocol": c}
    save(out / "stage-result.json", result)
    try:
        assert (
            sha(__file__) == spec["runtime_sha256"]
            and sha(out / "in.dynamics") == spec["input_sha256"]
        )
        for v in spec["sources"].values():
            assert sha(v["path"]) == v["sha256"]
        info = json.loads(Path(c["source_info"]).read_text())
        reference = np_geometry(c["source_data"], info)
        pre = out / "preflight"
        pre.mkdir()
        text = (out / "in.dynamics").read_text()
        (pre / "in.preflight").write_text(
            re.sub(r"^run \d+", "run 0", text, flags=re.M)
        )
        with (pre / "stdout.log").open("w") as f:
            subprocess.run(
                [
                    "srun",
                    "--ntasks=" + str(c["tasks"]),
                    c["lammps"],
                    "-in",
                    "in.preflight",
                    "-log",
                    "lammps.log",
                ],
                cwd=pre,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )
        with (out / "stdout.log").open("w") as f:
            subprocess.run(
                [
                    "srun",
                    "--ntasks=" + str(c["tasks"]),
                    c["lammps"],
                    "-in",
                    "in.dynamics",
                    "-log",
                    "lammps.log",
                ],
                cwd=out,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )
        log = (out / "lammps.log").read_text()
        if (
            re.search(r"^\s*ERROR(?:\s|:)|\b(?:nan|inf)\b", log, re.I | re.M)
            or "SILICONELAB_COMPOSITE_DYNAMICS_COMPLETE" not in log
        ):
            raise ValueError("LAMMPS error or missing completion")
        if (
            sum(map(int, re.findall(r"Loop time of .*? for (\d+) steps", log)))
            != c["dw_steps"] + c["msd_steps"]
        ):
            raise ValueError("Wrong integrated duration")
        result["trajectories"] = {
            k: frames(
                out / ("dump." + k + ".lammpstrj"),
                spec["expected_atoms"],
                (
                    range(0, c["dw_steps"] + 1, c["dw_stride"])
                    if k == "dw"
                    else range(
                        c["dw_steps"],
                        c["dw_steps"] + c["msd_steps"] + 1,
                        c["msd_stride"],
                    )
                ),
            )
            for k in ("dw", "msd")
        }
        final = np_geometry(out / "data.dynamics_final", info)
        error = max(
            abs(a - b) for m in reference for a, b in zip(reference[m], final[m])
        )
        if error > c["final_rigidity_tolerance_A"]:
            raise ValueError("Rigid NP geometry changed")
        for v in spec["sources"].values():
            assert sha(v["path"]) == v["sha256"]
        result.update(
            status="complete",
            maximum_final_NP_distance_error_A=error,
            final_data_sha256=sha(out / "data.dynamics_final"),
            warnings=sorted(set(s for s in log.splitlines() if "WARNING:" in s)),
        )
    except Exception as e:
        result.update(status="failed", error=str(e))
        raise
    finally:
        save(out / "stage-result.json", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["prepare", "run"])
    p.add_argument("--config", type=Path, required=True)
    a = p.parse_args()
    (prepare if a.action == "prepare" else run)(a.config)
