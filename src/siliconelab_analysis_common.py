"""Shared metadata/provenance and output contracts for parts 2 and 3."""

import csv
import hashlib
import importlib.machinery
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def static_driver():
    here = Path(__file__).resolve().parent
    p = next(
        (
            p
            for p in [here / "siliconelab_analyzer.py", here / "siliconelab_analyzer"]
            if p.exists()
        ),
        None,
    )
    if not p:
        raise RuntimeError("Install the part1 analyzer beside parts2/3")
    loader = importlib.machinery.SourceFileLoader("_siliconelab_static", str(p))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        return
    with Path(path).open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def backend(name, explicit=None):
    here = Path(__file__).resolve().parent
    for d in ([Path(explicit)] if explicit else [here, here.parent / "build"]):
        p = d / name
        if p.is_file():
            return str(p.resolve())
    raise RuntimeError(
        "Missing backend " + name + "; build/install C++ targets or set --backend-dir"
    )


def invoke(command, out):
    with (Path(out) / "commands.log").open("a") as f:
        f.write(json.dumps(command) + "\n")
        f.flush()
        subprocess.run(command, check=True, stdout=f, stderr=subprocess.STDOUT)


def metadata(info, geometry=None):
    source = json.loads(Path(info).read_text())
    return source, static_driver().normalize(source, geometry)


def state_row(name, status, reason="", files=None):
    return {"analysis": name, "status": status, "reason": reason, "files": files or []}


def report_markdown(out, title, rows, notes):
    text = (
        "# " + title + "\n\n| Analysis | Status | Evidence / reason |\n|---|---|---|\n"
    )
    for r in rows:
        text += (
            "| "
            + r["analysis"]
            + " | "
            + r["status"]
            + " | "
            + (r["reason"] or ", ".join(r["files"])).replace("|", "/")
            + " |\n"
        )
    text += "\n" + "\n\n".join(notes) + "\n"
    (Path(out) / "report.md").write_text(text)
