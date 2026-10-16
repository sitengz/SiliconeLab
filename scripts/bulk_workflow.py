#!/usr/bin/env python3
"""Prepare, execute and submit bulk formation/tensile chains with saved contracts.

Generation and validation belong on compute nodes. Submission is an explicit,
separate action. Tensile scheduling is serial; every axis reads formation data.
"""
import argparse
import base64
import collections
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zlib

ROOT = Path(__file__).resolve().parents[1]
LAMMPS_MODULE = "lammps/20230802.2-py310-openmpi4-ezoqd7f"
STAGES = ("formation", "tensile-x", "tensile-y", "tensile-z")


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def commands(text):
    return [line.split("#", 1)[0].split() for line in text.splitlines()
            if line.split("#", 1)[0].strip()]


def has_lammps_error(text):
    return bool(re.search(r"^\s*ERROR(?:\s|:)|\b(?:nan|inf)\b", text, re.I | re.M))


def cold_header(text, data):
    """Use the actual final native stage coefficients, including scaled NP pairs."""
    settings, bonded, pairs = {}, {}, {}
    for w in commands(text):
        if w[0] in ("units", "boundary", "atom_style", "bond_style", "angle_style",
                    "dihedral_style", "special_bonds", "comm_modify", "pair_style"):
            settings[w[0]] = w
            if w[0] == "pair_style":
                pairs = {}
        elif w[0] in ("bond_coeff", "angle_coeff", "dihedral_coeff"):
            bonded[(w[0], w[1])] = w
        elif w[0] == "pair_coeff":
            pairs[tuple(w[1:3])] = w
    if settings["boundary"] != ["boundary", "p", "p", "p"]:
        raise ValueError("Bulk tensile requires periodic bulk")
    if settings["pair_style"] != ["pair_style", "lj/gromacs", "12", "15"]:
        raise ValueError("This bulk workflow requires the native cold PDMS pair potential")
    order = ("units", "boundary", "atom_style", "bond_style", "angle_style", "dihedral_style",
             "special_bonds", "pair_style", "comm_modify")
    return "\n".join(" ".join(settings[k]) for k in order) + "\nread_data " + data + "\n" + \
        "\n".join(" ".join(w) for w in [*bonded.values(), *pairs.values()]) + "\n"


def tensile_input(text, source, axis, composite, rate, strain):
    steps_float = strain / rate / 5
    steps = round(steps_float)
    if not math.isclose(steps_float, steps, abs_tol=1e-6) or steps % 1000:
        raise ValueError("Loading duration must fit complete 1000-step sampling blocks")
    trans = [a for a in "xyz" if a != axis]
    index = "xyz".index(axis) + 1
    pressure = " ".join(a + " 1 1 500" for a in trans)
    area = "*".join("l" + a for a in trans)
    reference = "*".join("v_L" + a + "0" for a in trans)
    np_group = "group nanoparticles type 6\ngroup precursor subtract all nanoparticles\nneigh_modify exclude molecule/intra nanoparticles\n" if composite else "group precursor type 1 2 3 4 5\n"
    rigid = "fix np_integrate nanoparticles rigid/nvt/small molecule temp 300 300 50 reinit no\n" if composite else ""
    return f"""# Independent {axis} loading from the common formation snapshot; scheduling is serial.
{cold_header(text, source)}
{np_group}neighbor 2 bin
neigh_modify delay 5 every 1
compute precursor_temp precursor temp
compute composite_temp all temp
compute composite_press all pressure composite_temp
timestep 5
reset_timestep 0
thermo 1000
thermo_style custom step time c_precursor_temp c_composite_temp density lx ly lz c_composite_press[1] c_composite_press[2] c_composite_press[3] etotal epair ebond eangle edihed
thermo_modify lost error flush yes
{rigid}fix initial precursor nvt temp 300 300 50
run 0
unfix initial
variable Lx0 equal $(lx:%.17g)
variable Ly0 equal $(ly:%.17g)
variable Lz0 equal $(lz:%.17g)
variable erate equal {rate:.17g}
print "# Lx0_A Ly0_A Lz0_A" file reference_box.dat screen no
print "${{Lx0}} ${{Ly0}} ${{Lz0}}" append reference_box.dat screen no
variable t equal time
variable T equal c_precursor_temp
variable Tall equal c_composite_temp
variable ex equal lx/v_Lx0-1
variable ey equal ly/v_Ly0-1
variable ez equal lz/v_Lz0-1
variable px equal c_composite_press[1]
variable py equal c_composite_press[2]
variable pz equal c_composite_press[3]
variable Lx equal lx
variable Ly equal ly
variable Lz equal lz
variable V equal vol
variable true_stress equal -c_composite_press[{index}]*0.101325
variable engineering_stress equal v_true_stress*({area})/({reference})
fix integrate precursor npt temp 300 300 50 {pressure} couple none dilate precursor
fix_modify integrate temp precursor_temp press composite_press
fix stretch precursor deform 1 {axis} erate ${{erate}} units box remap x
fix record all ave/time 100 10 1000 v_t v_T v_Tall v_ex v_ey v_ez v_px v_py v_pz v_Lx v_Ly v_Lz v_V v_true_stress v_engineering_stress file stress_strain.dat
dump geom all custom 100000 dump.tensile.lammpstrj id mol type x y z ix iy iz
dump_modify geom first yes sort id format float %.12g
restart 500000 restart.tensile.1 restart.tensile.2
run {steps}
unfix record
undump geom
unfix stretch
unfix integrate
fix hold precursor nvt temp 300 300 50
run 1 start 0 stop $(step+1:%.0f)
write_data data.tensile_final nocoeff
write_restart restart.tensile_final
write_dump all custom dump.tensile_final id mol type x y z ix iy iz modify sort id format float %.12g
print "SILICONELAB_STAGE_COMPLETE tensile-{axis}" file completion.txt screen yes
unfix hold
{"unfix np_integrate" if composite else ""}
"""


def initialize(batch):
    batch = Path(batch).resolve()
    batch.mkdir(parents=True, exist_ok=False)
    cases = []
    for label, lam in [("pure-v22", None), *[("v22-np32-lambda-" + f"{i / 10:.1f}".replace(".", "p"), i / 10) for i in range(10, 4, -1)]]:
        case = batch / label
        (case / "configs").mkdir(parents=True)
        config = ("system = coating\nmodel = v22\n" if lam is None else "system = composite\nbase = v22\n")
        config += f"output_dir = ../generated/{label}\n"
        config += "n1 = 128\nm1 = 900\nn2 = 32\nfunctionality = 8\nn4 = 5\nm4 = 6\ndensity = 0.1\ntarget_density = 0.8\ncrosslink_distribution = random\ncrosslink_seed = 20260722\nseed = 5489\n"
        if lam is not None:
            config += f"np_count = 32\nnp_radius = 30\nnp_seed = 20261014\nnp_interaction = native_pdms\nnp_lambda = {lam:.1f}\n"
        config += "tensile_erate = 2e-8\ntensile_strain = 1\n"
        (case / "configs/model.conf").write_text(config)
        cases.append({"id": label, "lambda": lam, "np_count": 0 if lam is None else 32,
                      "expected_atoms": 122430 if lam is None else 136158})
    plan = {"schema_version": 1, "batch": batch.name, "cases": cases,
            "stages": list(STAGES), "source_rule": "Each axis starts from formation; x/y/z dependencies serialize scheduling only",
            "independent_job_limit": 8, "resources": {"tasks": 96, "memory": "200G", "time": "48:00:00"},
            "tensile": {"engineering_rate_fs_inverse": 2e-8, "strain": 1, "steps": 10000000,
                        "temperature_K": 300, "transverse_pressure_atm": 1, "output_consistency_steps": 1},
            "automatic_MD_resubmission": False}
    save(batch / "workflow.json", plan)
    save(batch / "workflow-state.json", {"status": "planned", "created_at": now(), "stages": {}})
    return plan


def launcher(directory, controller=None):
    controller = controller or ROOT / "scripts/bulk_workflow.py"
    return f"""#!/bin/bash
#SBATCH --account=wxia
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --chdir={directory}
#SBATCH --output={directory}/%j.out
#SBATCH --error={directory}/%j.err
set -eo pipefail
source /etc/profile
source /work/wxia/siteng/log/tools/nova-logging.sh
nova_run 'Load established simulation module' --summary 'Loaded validated LAMMPS module and its native MPI dependencies' module load {LAMMPS_MODULE}
export OMP_NUM_THREADS=1
nova_run 'Execute and verify bulk formation or tensile stage' --summary 'Full simulation output and verification retained in this stage folder; stage-result.json gates dependent work' python3 {controller} run --stage-dir {directory}
"""


def prepare(batch, generator):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Production generation requires a compute allocation")
    batch = Path(batch).resolve()
    plan = load(batch / "workflow.json")
    state = load(batch / "workflow-state.json")
    if state["status"] != "planned":
        raise ValueError("Preserve existing preparation and run evidence")
    state.update(status="preparing", preparation_job_id=os.environ["SLURM_JOB_ID"])
    save(batch / "workflow-state.json", state)
    (batch / "runtime").mkdir()
    controller = batch / "runtime/bulk_workflow.py"
    shutil.copy2(Path(__file__), controller)
    shared_atoms = None
    for item in plan["cases"]:
        case = batch / item["id"]
        config = case / "configs/model.conf"
        subprocess.run([str(generator), "--config", str(config)], check=True)
        package = case / "generated" / item["id"]
        manifest = load(package / "siliconelab.json")
        if item["lambda"] is not None:
            data = (package / manifest["data_file"]).read_text()
            invariant = data.split("\n", 1)[1]
            if shared_atoms is not None and invariant != shared_atoms:
                raise RuntimeError("Lambda cases must start from identical matrix/NP coordinates and topology")
            shared_atoms = invariant
        original = (package / ("in." + item["id"])).read_text()
        for index, stage in enumerate(STAGES):
            directory = case / stage
            directory.mkdir()
            scientific_parent = None if stage == "formation" else str(case / "formation")
            final = "data." + item["id"] + ".formation_final" if stage == "formation" else "data.tensile_final"
            if stage == "formation":
                for p in package.iterdir():
                    if p.is_file():
                        shutil.copy2(p, directory / p.name)
                inp = "in." + item["id"]
                text = original + f"\nwrite_data {final} nocoeff\nwrite_restart restart.formation_final\nprint \"SILICONELAB_STAGE_COMPLETE formation\" file completion.txt screen yes\n"
            else:
                inp = "in.tensile"
                text = tensile_input(original, "../formation/data." + item["id"] + ".formation_final",
                                     stage[-1], item["lambda"] is not None,
                                     manifest["bulk_tensile"]["engineering_rate_fs_inverse"],
                                     manifest["bulk_tensile"]["target_engineering_strain"])
            (directory / inp).write_text(text)
            (directory / "submit.sbatch").write_text(launcher(directory, controller))
            spec = {**item, "stage": stage, "input": inp, "input_sha256": sha(directory / inp),
                    "config": str(config), "config_sha256": sha(config), "package": str(package),
                    "initial_data": manifest["data_file"], "initial_data_sha256": sha(package / manifest["data_file"]),
                    "final_data": final, "scientific_parent": scientific_parent,
                    "executor_sha256": sha(controller),
                    "scheduling_parent": None if index == 0 else str(case / STAGES[index - 1]),
                    "expected_MD_steps": 8000000 + (3 if item["lambda"] is not None else 0) if stage == "formation" else 10000001,
                    "source_revision": subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()}
            save(directory / "stage.json", spec)
            state["stages"][item["id"] + ":" + stage] = {"directory": str(directory), "state": "PREPARED", "job_id": None}
            save(batch / "workflow-state.json", state)
    state.update(status="prepared_inputs", finished_preparing_at=now())
    save(batch / "workflow-state.json", state)


def snapshot(path, expected, np_start=None):
    """Read the complete atom section, including rigid-body unwrapped positions."""
    count, atoms, box, section = 0, {}, {}, None
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
                ident, mol, typ = map(int, w[:3])
                if ident in atoms or not 1 <= ident <= expected or len(w) not in (7, 10):
                    raise ValueError("Invalid final atom identity")
                xyz = [float(w[4 + j]) + (int(w[7 + j]) * box[a] if len(w) == 10 else 0) for j, a in enumerate("xyz")]
                if not all(map(math.isfinite, xyz)):
                    raise ValueError("Nonfinite final atom coordinate")
                atoms[ident] = (mol, typ, xyz)
    if count != expected or len(atoms) != expected or len(box) != 3:
        raise ValueError("Snapshot count/box contract failed")
    return atoms, box


def rigid_reference(spec):
    if spec["lambda"] is None:
        return {}
    atoms, _ = snapshot(Path(spec["package"]) / spec["initial_data"], spec["expected_atoms"])
    groups = collections.defaultdict(list)
    for ident, (mol, typ, xyz) in atoms.items():
        if typ == 6:
            groups[mol].append((ident, xyz))
    if len(groups) != 32 or any(len(rows) != 429 for rows in groups.values()):
        raise ValueError("NP body identity/count contract failed")
    return {mol: {"ids": [r[0] for r in rows], "anchors": [r[0] for r in rows[:4]],
                  "distances": [[math.dist(a[1], r[1]) for r in rows] for a in rows[:4]]}
            for mol, rows in groups.items()}


def rigidity_error(positions, reference):
    worst = 0.0
    for body in reference.values():
        for anchor, distances in zip(body["anchors"], body["distances"]):
            a = positions[anchor]
            for ident, distance in zip(body["ids"], distances):
                worst = max(worst, abs(math.dist(a, positions[ident]) - distance))
    return worst


def scan_trajectory(path, expected_atoms, expected_steps, reference):
    steps, worst = [], 0.0
    selected = {i for b in reference.values() for i in b["ids"]}
    with Path(path).open() as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError("Incomplete trajectory header")
            step = int(f.readline())
            if f.readline().strip() != "ITEM: NUMBER OF ATOMS" or int(f.readline()) != expected_atoms:
                raise ValueError("Trajectory atom count mismatch")
            if f.readline().strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError("Unexpected trajectory boundary")
            lengths = [float((w := f.readline().split())[1]) - float(w[0]) for _ in range(3)]
            columns = f.readline().split()[2:]
            ix = [columns.index(a) for a in ("id", "type", "x", "y", "z", "ix", "iy", "iz")]
            ids, positions = set(), {}
            for _ in range(expected_atoms):
                w = f.readline().split()
                if len(w) != len(columns):
                    raise ValueError("Incomplete atom record in trajectory")
                ident = int(w[ix[0]])
                if ident in ids or not 1 <= ident <= expected_atoms:
                    raise ValueError("Duplicate/missing trajectory atom IDs")
                ids.add(ident)
                xyz = [float(w[ix[2 + j]]) for j in range(3)]
                if not all(map(math.isfinite, xyz)):
                    raise ValueError("Nonfinite trajectory coordinates")
                if ident in selected:
                    positions[ident] = [xyz[j] + int(w[ix[5 + j]]) * lengths[j] for j in range(3)]
            if reference:
                worst = max(worst, rigidity_error(positions, reference))
            steps.append(step)
    if steps != list(expected_steps) or worst > .02:
        raise ValueError("Trajectory duration/sampling or NP rigidity contract failed")
    return {"frames": len(steps), "first_step": steps[0], "last_step": steps[-1], "maximum_NP_distance_error_A": worst}


def validate_stage(directory, spec, result):
    atoms, box = snapshot(directory / spec["final_data"], spec["expected_atoms"])
    reference = rigid_reference(spec)
    error = rigidity_error({i: r[2] for i, r in atoms.items() if r[1] == 6}, reference)
    if error > 1e-6:
        raise ValueError("Final NP geometry changed")
    text = (directory / "lammps.log").read_text()
    if has_lammps_error(text) or ("SILICONELAB_STAGE_COMPLETE " + spec["stage"]) not in text:
        raise ValueError("LAMMPS error/nonfinite output or missing completion marker")
    steps = sum(map(int, re.findall(r"Loop time of .*? for (\d+) steps", text)))
    if steps != spec["expected_MD_steps"]:
        raise ValueError(f"Executed {steps} steps, expected {spec['expected_MD_steps']}")
    if spec["stage"] == "formation":
        trajectory = scan_trajectory(directory / "dump.msd.lammpstrj", spec["expected_atoms"], range(0, 1000001, 1000), reference)
    else:
        trajectory = scan_trajectory(directory / "dump.tensile.lammpstrj", spec["expected_atoms"], range(0, 10000001, 100000), reference)
        rows = [list(map(float, line.split())) for line in (directory / "stress_strain.dat").read_text().splitlines() if line.strip() and not line.startswith("#")]
        axis_column = 4 + "xyz".index(spec["stage"][-1])
        if len(rows) != 10000 or any(len(r) != 16 or int(r[0]) != (i + 1) * 1000 or not all(map(math.isfinite, r)) for i, r in enumerate(rows)):
            raise ValueError("Incomplete stress sampling")
        if max(abs(r[axis_column] - (r[0] * 5 - 2250) * 2e-8) for r in rows) > 1e-5:
            raise ValueError("Tensile strain/clock contract failed")
        result["stress_samples"] = len(rows)
    result.update(final_atoms=len(atoms), final_box_lengths_A=box, final_data_sha256=sha(directory / spec["final_data"]),
                  executed_steps=steps, trajectory=trajectory, maximum_final_NP_distance_error_A=error,
                  warnings=sorted(set(line.strip() for line in text.splitlines() if "WARNING:" in line)))


def run_stage(directory):
    directory = Path(directory).resolve()
    spec = load(directory / "stage.json")
    result_path = directory / "stage-result.json"
    if result_path.exists() or (directory / "lammps.log").exists():
        raise ValueError("Existing stage evidence must be preserved; no automatic resubmission")
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    result = {"status": "running", "stage": spec["stage"], "case": spec["id"], "job_id": os.environ["SLURM_JOB_ID"], "started_at": now()}
    save(result_path, result)
    try:
        if sha(Path(__file__)) != spec["executor_sha256"]:
            raise ValueError("Pinned stage executor changed")
        if sha(directory / spec["input"]) != spec["input_sha256"] or sha(spec["config"]) != spec["config_sha256"]:
            raise ValueError("Input/model configuration changed since preparation")
        if sha(Path(spec["package"]) / spec["initial_data"]) != spec["initial_data_sha256"]:
            raise ValueError("Generated source changed")
        for parent in {spec["scientific_parent"], spec["scheduling_parent"]} - {None}:
            contract = load(Path(parent) / "stage-result.json")
            if contract["status"] != "complete":
                raise ValueError("Parent contract incomplete")
            parent_spec = load(Path(parent) / "stage.json")
            if sha(Path(parent) / parent_spec["final_data"]) != contract["final_data_sha256"]:
                raise ValueError("Parent data hash mismatch")
        if spec["scientific_parent"]:
            contract = load(Path(spec["scientific_parent"]) / "stage-result.json")
            result["formation_source_sha256"] = contract["final_data_sha256"]
        subprocess.run(["srun", "--ntasks=96", "lmp", "-in", spec["input"], "-log", "lammps.log"], cwd=directory, check=True)
        validate_stage(directory, spec, result)
        result.update(status="complete", finished_at=now())
    except Exception as e:
        result.update(status="failed", error=str(e), finished_at=now())
        raise
    finally:
        save(result_path, result)


def submit(batch):
    batch = Path(batch).resolve()
    state = load(batch / "workflow-state.json")
    if state["status"] not in ("validated", "submitted"):
        raise ValueError("Compute-node input preflight/smoke validation is required")
    for item in load(batch / "workflow.json")["cases"]:
        predecessor = None
        for stage in STAGES:
            key = item["id"] + ":" + stage
            entry = state["stages"][key]
            if entry.get("job_id"):
                predecessor = str(entry["job_id"])
                continue
            queue = subprocess.check_output(["squeue", "--array", "-u", os.environ["USER"], "-h", "-o", "%i|%E|%j"], text=True).splitlines()
            roots = [r for r in queue if r.split("|")[1] in ("", "(null)", "N/A")]
            if predecessor is None and len(roots) >= 8:
                raise RuntimeError("Eight independent queue elements reached; saved accepted jobs, resume later")
            name = "SLbulk-" + item["id"] + "-" + stage
            if any(r.split("|")[-1] == name for r in queue):
                raise RuntimeError("Matching unrecorded job exists; reconcile before retrying submission")
            history = subprocess.check_output(["sacct", "-S", state["created_at"][:10], "--name", name, "-n", "-P", "--format=JobIDRaw"], text=True).strip()
            if history:
                raise RuntimeError("Historical matching job exists; reconcile before submission")
            args = ["sbatch", "--parsable", "--job-name=" + name, "--chdir=" + entry["directory"]]
            if predecessor:
                args += ["--dependency=afterok:" + predecessor, "--kill-on-invalid-dep=yes"]
            args += [str(Path(entry["directory"]) / "submit.sbatch")]
            accepted = subprocess.check_output(args, text=True).strip().split(";", 1)[0]
            if not accepted.isdigit():
                raise RuntimeError("Unparseable accepted job ID; inspect scheduler before retrying")
            entry.update(job_id=accepted, state="SUBMITTED", dependency_job_id=predecessor, submitted_at=now(), submission_command=args)
            save(batch / "workflow-state.json", state)
            predecessor = accepted
            print("SUMMARY: Accepted " + key + " as " + accepted, flush=True)
    state.update(status="submitted", submitted_at=now())
    save(batch / "workflow-state.json", state)


def check(batch):
    batch = Path(batch).resolve()
    state = load(batch / "workflow-state.json")
    ids = [v["job_id"] for v in state["stages"].values() if v.get("job_id")]
    if not ids:
        raise ValueError("No submitted jobs to monitor")
    accounting = subprocess.check_output(["sacct", "-n", "-P", "-j", ",".join(ids),
                                         "--format=JobIDRaw,State,ExitCode,Elapsed,Timelimit"], text=True)
    records = {r[0]: r for line in accounting.splitlines() if (r := line.split("|")) and r[0] in ids}
    queue = subprocess.check_output(["squeue", "--array", "-u", os.environ["USER"], "-h", "-o", "%i|%T|%E|%j"], text=True)
    running = sum(line.split("|")[1] == "RUNNING" for line in queue.splitlines())
    if running > 8:
        raise RuntimeError("Account has more than eight running elements; inspect unrelated submissions")
    counts = collections.Counter()
    for key, entry in state["stages"].items():
        row = records.get(entry.get("job_id"))
        if not row:
            raise RuntimeError("Accounting record missing: " + key)
        entry.update(state=row[1], exit_code=row[2], elapsed=row[3], time_limit=row[4], last_checked_at=now())
        if row[4] not in ("2-00:00:00", "48:00:00"):
            raise RuntimeError("Unexpected production time limit: " + key)
        contract = Path(entry["directory"]) / "stage-result.json"
        if contract.exists():
            entry["result"] = load(contract)
        if row[1] == "COMPLETED" and (row[2] != "0:0" or entry.get("result", {}).get("status") != "complete"):
            raise RuntimeError("Completed scheduler job lacks a valid stage contract: " + key)
        counts[row[1]] += 1
    state.update(last_checked_at=now(), account_running_elements=running,
                 all_terminal=all(v["state"].split()[0] in ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL") for v in state["stages"].values()))
    save(batch / "workflow-state.json", state)
    print("SUMMARY: Bulk interface series: " + ", ".join(str(v) + " " + k.lower() for k, v in sorted(counts.items())) + ".")
    print("SUMMARY: Account running elements: " + str(running) + "; maximum eight; scientific stage results remain in each case folder.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("init", "prepare", "run", "submit", "check", "snapshot"))
    p.add_argument("--batch-dir", type=Path)
    p.add_argument("--stage-dir", type=Path)
    p.add_argument("--generator", type=Path, default=ROOT / "build/siliconelab_generator")
    a = p.parse_args()
    if a.action == "init":
        initialize(a.batch_dir)
    elif a.action == "prepare":
        prepare(a.batch_dir, a.generator)
    elif a.action == "run":
        run_stage(a.stage_dir)
    elif a.action == "submit":
        submit(a.batch_dir)
    elif a.action == "check":
        check(a.batch_dir)
    elif a.action == "snapshot":
        content = {"workflow": load(a.batch_dir / "workflow.json"), "state": load(a.batch_dir / "workflow-state.json"),
                   "preparation": load(a.batch_dir / "preparation-result.json")}
        print("SNAPSHOT:" + base64.b64encode(zlib.compress(json.dumps(content).encode())).decode())


if __name__ == "__main__":
    main()
