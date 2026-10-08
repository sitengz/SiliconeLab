#!/usr/bin/env python3
"""Submit exactly two metadata-correlated analysis jobs for each prepared film."""

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import film_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--installation", type=Path, required=True)
    p.add_argument("--python-lib", type=Path, required=True)
    p.add_argument("--validation-job", required=True)
    a = p.parse_args()
    state = w.load(w.STATE)
    if not state.get("preparation_complete"):
        raise RuntimeError("Film preparation incomplete")
    state.setdefault("analysis_jobs", {})
    for case in state["cases"]:
        for part in ["part2", "part3"]:
            key = case + ":" + part
            if state["analysis_jobs"].get(key, {}).get("job_id"):
                continue
            parents = [a.validation_job, state["stages"][case + ":dynamics"]["job_id"]]
            if part == "part3":
                parents.append(state["analysis_jobs"][case + ":part2"]["job_id"])
                if state["cases"][case]["system"] != "oil":
                    parents += [
                        state["stages"][case + ":tensile-" + axis]["job_id"]
                        for axis in "xy"
                    ]
            directory = w.ROOT / "runs" / case / ("analysis-job-" + part)
            directory.mkdir(exist_ok=True)
            output = w.ROOT / "runs" / case / ("analysis-" + part)
            if output.exists():
                raise RuntimeError(
                    "Existing analysis results; refuse overwrite " + str(output)
                )
            name = "SLfilm-" + case + "-" + part
            queue = subprocess.check_output(
                ["squeue", "--array", "-u", os.environ["USER"], "-h", "-o", "%i|%E|%j"],
                text=True,
            ).splitlines()
            if any(r.split("|")[-1] == name for r in queue):
                raise RuntimeError("Matching unrecorded queued analysis: " + name)
            history = subprocess.check_output(
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
            if history:
                raise RuntimeError(
                    "Matching historical analysis requires reconciliation: " + name
                )
            command = [
                "python3",
                str(w.ROOT / "scripts/analyze_film.py"),
                part,
                case,
                "--installation",
                str(a.installation),
                "--output-dir",
                str(output),
            ]
            script = (
                "#!/bin/bash\n#SBATCH --account=wxia\n#SBATCH --partition=nova\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=2\n#SBATCH --mem=32G\n#SBATCH --time=08:00:00\nset -eo pipefail\nsource /etc/profile\nsource /work/wxia/siteng/log/tools/nova-logging.sh\nmodule load cmake/3.31.8-sds7j4t\nexport OMP_NUM_THREADS=2\nexport PYTHONPATH="
                + shlex.quote(str(a.python_lib))
                + "${PYTHONPATH:+:$PYTHONPATH}\nnova_run "
                + shlex.quote("Analyze " + case + " " + part)
                + " --summary "
                + shlex.quote(
                    "Write correlated scientific report with explicit source, applicability and quality statuses; full output remains in the analysis job folder."
                )
                + " "
                + shlex.join(command)
                + "\n"
            )
            (directory / "analyze.sbatch").write_text(script)
            cmd = [
                "sbatch",
                "--parsable",
                "--job-name=" + name,
                "--chdir=" + str(directory),
                "--output=" + str(directory / "%j.out"),
                "--error=" + str(directory / "%j.err"),
                "--dependency=afterok:" + ":".join(parents),
                "--kill-on-invalid-dep=yes",
                "analyze.sbatch",
            ]
            job = (
                subprocess.check_output(cmd, cwd=directory, text=True)
                .strip()
                .split(";")[0]
            )
            if not job.isdigit():
                raise RuntimeError("Unexpected sbatch reply")
            item = {
                "case": case,
                "part": part,
                "job_id": job,
                "state": "SUBMITTED",
                "directory": str(directory.relative_to(w.ROOT)),
                "output_directory": str(output.relative_to(w.ROOT)),
                "dependency_job_ids": parents,
                "installation": str(a.installation),
                "python_lib": str(a.python_lib),
                "slurm_stdout": str(directory / (job + ".out")),
                "slurm_stderr": str(directory / (job + ".err")),
                "submitted_at": w.now(),
                "submission_command": cmd,
            }
            state["analysis_jobs"][key] = item
            w.save(directory / "submission.json", item)
            w.save(w.STATE, state)
            print(
                "SUMMARY: "
                + key
                + " job "
                + job
                + " accepted; dependency jobs "
                + ",".join(parents),
                flush=True,
            )


if __name__ == "__main__":
    main()
