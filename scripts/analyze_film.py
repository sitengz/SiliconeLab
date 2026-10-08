#!/usr/bin/env python3
"""Compute-node analysis runner gated by the film workflow's stage contracts."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import film_workflow as workflow

ROOT = workflow.ROOT


def main():
    p = argparse.ArgumentParser()
    p.add_argument("part", choices=["part2", "part3"])
    p.add_argument("case")
    p.add_argument("--installation", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Analysis must run on compute nodes")
    state = workflow.load(workflow.STATE)
    c = state["cases"][a.case]
    run = ROOT / "runs" / a.case
    info = ROOT / c["source_info"]
    if a.part == "part2":
        item = state["stages"][a.case + ":dynamics"]
        workflow.check_stage(item)
        directory = ROOT / item["directory"]
        command = [
            sys.executable,
            str(a.installation / "bin/siliconelab_part2"),
            str(directory / "data.dynamics_final"),
            str(info),
            "--geometry",
            "film",
            "--trajectory",
            str(directory / "dump.layer_dynamics.lammpstrj"),
            "--dw-trajectory",
            str(directory / "dump.debye_waller.lammpstrj"),
            "--output-dir",
            str(a.output_dir),
        ]
    else:
        stage = "film" if c["system"] == "oil" else "surface"
        item = state["stages"][a.case + ":" + stage]
        workflow.check_stage(item)
        directory = ROOT / item["directory"]
        prefix = (
            "energy." + a.case + ".film"
            if c["system"] == "oil"
            else "pressure.production"
        )
        prod = directory / (prefix + ".dat")
        equil = directory / (
            "energy." + a.case + ".film_eq.dat"
            if c["system"] == "oil"
            else "pressure.equil.dat"
        )
        command = [
            sys.executable,
            str(a.installation / "bin/siliconelab_part3"),
            str(info),
            "--geometry",
            "film",
            "--pressure-production",
            str(prod),
            "--pressure-equil",
            str(equil),
            "--initial-data",
            str(run / ("data." + a.case)),
            "--geometry-data",
            str(directory / item["final_data"]),
            "--part2-report",
            str(run / "analysis-part2" / "analysis.json"),
            "--output-dir",
            str(a.output_dir),
        ]
        if c["system"] != "oil":
            for axis in "xy":
                workflow.check_stage(state["stages"][a.case + ":tensile-" + axis])
                command += [
                    "--tensile-" + axis,
                    str(run / ("tensile-" + axis) / "stress_strain.dat"),
                ]
            command += ["--material-thickness-A", str(c["nominal_thickness_A"])]
    subprocess.run(command, cwd=a.output_dir.parent, check=True)
    report = json.loads((a.output_dir / "analysis.json").read_text())
    if report["status"] != "complete":
        raise RuntimeError("Analyzer failed")
    print(
        "SUMMARY: "
        + a.case
        + " "
        + a.part
        + " report written; per-observable applicability and quality statuses remain explicit."
    )


if __name__ == "__main__":
    main()
