#!/usr/bin/env python3
"""Film mechanical surface statistics, tensile response and oil ATSC4i."""

import argparse
import collections
import importlib.machinery
import importlib.util
import json
import math
from pathlib import Path
import statistics
from siliconelab_analysis_common import *

PRESSURE_FIELDS = [
    "time_fs",
    "temp_K",
    "pe_kcal_per_mol",
    "pxx_atm",
    "pyy_atm",
    "pzz_atm",
    "lx_A",
    "ly_A",
    "lz_A",
    "wall_lo_force",
    "wall_hi_force",
]
ANALYSES = [
    "oil_film_surface_tension",
    "surface_statistics_and_plots",
    "network_film_surface_stress",
    "tensile_response",
    "ATSC4i",
]


def numeric_rows(path, columns):
    result = []
    with Path(path).open() as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            words = line.split()
            if len(words) != len(columns):
                raise ValueError(
                    f"{path}:{n}: expected {len(columns)} columns, got {len(words)}"
                )
            values = list(map(float, words))
            if not all(map(math.isfinite, values)):
                raise ValueError("Nonfinite input table value")
            result.append(dict(zip(columns, values)))
    if not result:
        raise ValueError("Empty measurement table: " + str(path))
    return result


def regression(x, y):
    if len(x) < 3:
        return None
    mx = statistics.mean(x)
    my = statistics.mean(y)
    sx = sum((a - mx) ** 2 for a in x)
    if not sx:
        return None
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sx
    intercept = my - slope * mx
    rss = sum((b - intercept - slope * a) ** 2 for a, b in zip(x, y))
    sst = sum((b - my) ** 2 for b in y)
    return {
        "slope": slope,
        "intercept": intercept,
        "R2": 1 - rss / sst if sst else None,
        "points": len(x),
        "x_min": min(x),
        "x_max": max(x),
    }


def surface_table(path, out, block_ns=5, expected_ns=50, phase="production"):
    raw = numeric_rows(path, PRESSURE_FIELDS)
    times = [r["time_fs"] for r in raw]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError(
            "Non-increasing phase clock; do not concatenate equilibration and production"
        )
    if any(r["lz_A"] <= 0 for r in raw):
        raise ValueError("Invalid cell height")
    rows = []
    for r in raw:
        r = dict(r)
        r["time_ns"] = r["time_fs"] * 1e-6
        r["surface_mechanical_mN_m"] = (
            0.0101325
            * 0.5
            * r["lz_A"]
            * (r["pzz_atm"] - 0.5 * (r["pxx_atm"] + r["pyy_atm"]))
        )
        r["wall_contact"] = (
            abs(r["wall_lo_force"]) > 1e-10 or abs(r["wall_hi_force"]) > 1e-10
        )
        rows.append(r)
    write_tsv(out / ("surface_" + phase + "_series.tsv"), rows)
    # Samples occur at block endpoints; each (b*width,(b+1)*width] is disjoint.
    grouped = collections.defaultdict(list)
    for r in rows:
        grouped[max(0, math.ceil(r["time_ns"] / block_ns - 1e-12) - 1)].append(r)
    blocks = []
    finished = max(r["time_ns"] for r in rows)
    for b, data in sorted(grouped.items()):
        complete = finished + 1e-9 >= (b + 1) * block_ns
        blocks.append(
            {
                "block": b,
                "start_ns": b * block_ns,
                "end_ns": (b + 1) * block_ns,
                "samples": len(data),
                "complete": complete,
                "mean_mN_m": statistics.mean(
                    r["surface_mechanical_mN_m"] for r in data
                ),
                "wall_contact_samples": sum(r["wall_contact"] for r in data),
            }
        )
    write_tsv(out / ("surface_" + phase + "_blocks.tsv"), blocks)
    means = [b["mean_mN_m"] for b in blocks if b["complete"]]
    sem = statistics.stdev(means) / math.sqrt(len(means)) if len(means) > 1 else None
    trend = regression(
        [0.5 * (b["start_ns"] + b["end_ns"]) for b in blocks if b["complete"]], means
    )
    result = {
        "status": (
            "invalid_wall_contact"
            if any(r["wall_contact"] for r in rows)
            else (
                "complete"
                if finished + 1e-9 >= expected_ns and len(means) >= 2
                else "insufficient_sampling"
            )
        ),
        "phase": phase,
        "duration_ns": finished,
        "expected_duration_ns": expected_ns,
        "samples": len(rows),
        "completed_blocks": len(means),
        "block_ns": block_ns,
        "sample_mean_mN_m": statistics.mean(r["surface_mechanical_mN_m"] for r in rows),
        "block_mean_mN_m": statistics.mean(means) if means else None,
        "block_SEM_mN_m": sem,
        "approx_normal_95pct_halfwidth_mN_m": 1.96 * sem if sem is not None else None,
        "wall_contact_samples": sum(r["wall_contact"] for r in rows),
        "block_drift_fit": trend,
        "method": "Lz(cell including vacuum)/2 * [Pzz-(Pxx+Pyy)/2]; two surfaces",
        "notes": [
            "Block SEM and normal confidence interval are approximate; choose blocks longer than correlations and inspect drift.",
            "No equilibration conclusion is inferred from elapsed time or a complete file.",
        ],
    }
    write_json(out / ("surface_" + phase + "_summary.json"), result)
    return rows, blocks, result


def tensile_table(path, axis, thickness, out):
    fields = [
        "step",
        "time_fs",
        "temp_K",
        "strain_x",
        "strain_y",
        "strain_z",
        "pxx_atm",
        "pyy_atm",
        "pzz_atm",
        "lx_A",
        "ly_A",
        "lz_A",
        "cell_volume_A3",
    ]
    raw = numeric_rows(path, fields)
    if thickness <= 0:
        raise ValueError("Positive material thickness required")
    if any(b["time_fs"] <= a["time_fs"] for a, b in zip(raw, raw[1:])):
        raise ValueError("Tensile time must increase")
    transverse = "y" if axis == "x" else "x"
    rows = []
    for r in raw:
        r = dict(r)
        material_volume = r["lx_A"] * r["ly_A"] * thickness
        correction = r["cell_volume_A3"] / material_volume
        for d in "xyz":
            r["sigma_" + d + "_MPa"] = -0.101325 * r["p" + d + d + "_atm"] * correction
        r["stress_difference_MPa"] = (
            r["sigma_" + axis + "_MPa"] - r["sigma_" + transverse + "_MPa"]
        )
        r["reduced_strain"] = (1 + r["strain_" + axis]) ** 2 - (
            1 + r["strain_" + transverse]
        ) ** 2
        rows.append(r)
    write_tsv(out / ("tensile_" + axis + "_series.tsv"), rows)
    fitrows = [r for r in rows if 0 <= r["strain_" + axis] <= 0.05]
    linear = regression(
        [r["strain_" + axis] for r in fitrows],
        [r["stress_difference_MPa"] for r in fitrows],
    )
    reduced = regression(
        [r["reduced_strain"] for r in fitrows],
        [r["stress_difference_MPa"] for r in fitrows],
    )
    end = max(r["strain_" + axis] for r in rows)
    result = {
        "status": "complete" if end >= 0.199 else "insufficient_loading",
        "axis": axis,
        "maximum_engineering_strain": end,
        "material_thickness_A": thickness,
        "volume_convention": "Lx*Ly*fixed nominal material thickness; vacuum corrected; actual thickness changes require a measured replacement",
        "small_strain_stress_difference_slope_MPa": linear,
        "neo_Hookean_reduced_strain_fit_MPa": reduced,
        "notes": [
            "Stress-difference slope is not automatically Young modulus.",
            "Reduced-strain fit is model dependent; nominal film thickness is an approximation.",
            "Two loading directions are independent samples of directional response, not independent replicas.",
        ],
    }
    write_json(out / ("tensile_" + axis + "_summary.json"), result)
    return rows, result


def descriptors(initial, normalized, out, descriptor_tool=None):
    candidates = (
        [Path(descriptor_tool)]
        if descriptor_tool
        else [
            Path(__file__).resolve().parent / "siliconelab_atsc4i.py",
            Path(__file__).resolve().parent.parent
            / "vendor/Silicone_Oil/Analysis/atsc4i.py",
        ]
    )
    p = next((p for p in candidates if p.exists()), None)
    if p is None:
        return {"status": "missing_dependency", "reason": "ATSC4i helper missing"}
    loader = importlib.machinery.SourceFileLoader("_atsc4i", str(p))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    # Scan the original generator IDs once; pendant type 5 is not a repeat.
    filler = normalized["components"]["filler"]
    sequences = collections.defaultdict(list)
    inside = False
    with initial.open() as f:
        for line in f:
            w = line.split()
            if w and w[0] == "Atoms":
                inside = True
                continue
            if inside and w:
                if not w[0].lstrip("-").isdigit():
                    break
                atom, mol, typ = map(int, w[:3])
                if filler["molecule_id_start"] <= mol <= filler[
                    "molecule_id_end"
                ] and typ in (1, 4):
                    sequences[mol].append((atom, "D" if typ == 1 else "M"))
    if len(sequences) != filler["M"]:
        raise ValueError("Original data does not cover all expected oil molecules")
    groups = collections.Counter(
        "".join(v for _, v in sorted(seq)) for seq in sequences.values()
    )
    values = []
    try:
        for sequence, count in sorted(groups.items()):
            values.append(dict(m.calculate_atsc4i(sequence), chain_count=count))
    except (ImportError, RuntimeError) as e:
        return {
            "status": "missing_dependency",
            "reason": str(e),
            "unique_sequences": len(groups),
            "oil_chains": sum(groups.values()),
            "terminal_model": "Provisional Me-[Si(R)(Rprime)-O]n-Si(Me)3; agree chemical end groups before comparison",
        }
    result = {
        "status": "complete",
        "oil_chains": sum(groups.values()),
        "unique_sequences": len(groups),
        "sequence_results": values,
        "terminal_model": "Provisional native capped graph; exact end groups and atomic weights must be agreed before comparison",
    }
    write_json(out / "ATSC4i.json", result)
    return result


def plot_surface(rows, blocks, path):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(
        [r["time_ns"] for r in rows],
        [r["surface_mechanical_mN_m"] for r in rows],
        alpha=0.3,
        lw=0.5,
        label="Samples",
    )
    ax.plot(
        [0.5 * (b["start_ns"] + b["end_ns"]) for b in blocks],
        [b["mean_mN_m"] for b in blocks],
        "o-",
        label="Block means",
    )
    ax.set(xlabel="Production time (ns)", ylabel="Mechanical surface estimate (mN/m)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return True


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("info", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--geometry", choices=["bulk", "film"])
    p.add_argument("--pressure-production", type=Path)
    p.add_argument("--pressure-equil", type=Path)
    p.add_argument("--block-ns", type=float, default=5)
    p.add_argument("--expected-production-ns", type=float, default=50)
    p.add_argument("--tensile-x", type=Path)
    p.add_argument("--tensile-y", type=Path)
    p.add_argument("--material-thickness-A", type=float)
    p.add_argument("--initial-data", type=Path)
    p.add_argument("--geometry-data", type=Path)
    p.add_argument("--backend-dir")
    p.add_argument("--part2-report", type=Path)
    p.add_argument("--descriptor-tool")
    p.add_argument("--no-plots", action="store_true")
    a = p.parse_args(argv)
    if a.output_dir.exists():
        raise ValueError("Choose a new output directory; results are never overwritten")
    if a.block_ns <= 0 or a.expected_production_ns <= 0:
        raise ValueError("Positive time windows required")
    source, normalized = metadata(a.info, a.geometry)
    out = a.output_dir
    out.mkdir(parents=True)
    rows = []
    result = {
        "status": "running",
        "source_info": str(a.info.resolve()),
        "source_info_sha256": sha(a.info),
        "analyses": rows,
        "source_files": {},
    }
    write_json(out / "analysis.json", result)
    try:
        oil = normalized["source_system"] == "oil"
        film = normalized["geometry"] == "film"
        surface = None
        geometry_check = {
            "status": "missing_output",
            "reason": "Surface geometry snapshot not provided",
        }
        if a.geometry_data and a.geometry_data.exists():
            import numpy as np
            from siliconelab_part2 import surfaces

            folder = out / "geometry"
            folder.mkdir()
            write_json(folder / "normalized.info", normalized)
            invoke(
                [
                    backend("siliconelab_spatial_dynamics", a.backend_dir),
                    "spatial",
                    str(a.geometry_data.resolve()),
                    str((folder / "normalized.info").resolve()),
                    str(folder.resolve()),
                ],
                out,
            )
            geometry_check = surfaces(
                read_tsv(folder / "component_profiles.tsv"), normalized["geometry"], np
            )
            result["source_files"][str(a.geometry_data)] = sha(a.geometry_data)
        result["geometry_check"] = geometry_check
        if a.part2_report:
            part2 = json.loads(a.part2_report.read_text())
            if part2.get("status") != "complete" or part2["source"][
                "info_sha256"
            ] != sha(a.info):
                raise RuntimeError(
                    "Part2 report must be complete and use the same source metadata"
                )
            result["correlated_part2_report"] = {
                "path": str(a.part2_report),
                "sha256": sha(a.part2_report),
                "surface": part2.get("surface"),
            }
        if a.pressure_production and a.pressure_production.exists() and film:
            series, blocks, surface = surface_table(
                a.pressure_production, out, a.block_ns, a.expected_production_ns
            )
            result["source_files"][str(a.pressure_production)] = sha(
                a.pressure_production
            )
            if a.pressure_equil and a.pressure_equil.exists():
                surface_table(
                    a.pressure_equil, out, a.block_ns, a.expected_production_ns, "equil"
                )
                result["source_files"][str(a.pressure_equil)] = sha(a.pressure_equil)
            plotted = (
                False
                if a.no_plots
                else plot_surface(series, blocks, out / "surface_production.png")
            )
            if (
                surface["status"] == "complete"
                and geometry_check["status"] != "complete"
            ):
                surface["status"] = "insufficient_geometry"
            result["surface"] = surface
            result["plot_status"] = (
                "disabled"
                if a.no_plots
                else ("complete" if plotted else "missing_matplotlib")
            )
        for i, applicable in [(0, oil and film), (2, not oil and film)]:
            rows.append(
                state_row(
                    ANALYSES[i],
                    (
                        surface["status"]
                        if applicable and surface
                        else ("missing_output" if applicable else "not_applicable")
                    ),
                    "Wall contact, block drift and film geometry must be checked",
                    (
                        ["surface_production_summary.json"]
                        if surface and applicable
                        else []
                    ),
                )
            )
        rows.append(
            state_row(
                ANALYSES[1],
                (
                    surface["status"]
                    if surface
                    else ("missing_output" if film else "not_applicable")
                ),
                "Statistics complete when sufficient blocks exist; plot availability is recorded separately",
                ["surface_production_blocks.tsv"] if surface else [],
            )
        )
        tensile = {}
        if not oil:
            thickness = a.material_thickness_A or normalized.get(
                "film_thickness_angstrom"
            )
            for axis, path in [("x", a.tensile_x), ("y", a.tensile_y)]:
                if path and path.exists():
                    if thickness is None:
                        raise ValueError(
                            "Tensile analysis requires --material-thickness-A or nominal film thickness metadata"
                        )
                    _, tensile[axis] = tensile_table(path, axis, float(thickness), out)
                    result["source_files"][str(path)] = sha(path)
        rows.append(
            state_row(
                ANALYSES[3],
                (
                    "not_applicable"
                    if oil
                    else (
                        "complete"
                        if len(tensile) == 2
                        and all(v["status"] == "complete" for v in tensile.values())
                        else ("partial" if tensile else "missing_output")
                    )
                ),
                "Independent x/y loading and a stated material-volume convention required",
                list("tensile_" + axis + "_summary.json" for axis in tensile),
            )
        )
        result["tensile"] = tensile
        hasoil = normalized["components"]["filler"]["M"] > 0
        if hasoil and a.initial_data and a.initial_data.exists():
            descriptor = descriptors(a.initial_data, normalized, out, a.descriptor_tool)
            result["source_files"][str(a.initial_data)] = sha(a.initial_data)
        else:
            descriptor = {
                "status": "missing_output" if hasoil else "not_applicable",
                "reason": (
                    "Original generated oil sequence required"
                    if hasoil
                    else "No oil chains"
                ),
            }
        rows.append(
            state_row(
                ANALYSES[4],
                descriptor["status"],
                descriptor.get(
                    "reason",
                    "Native atomistic reconstruction; end groups remain a stated chemical assumption",
                ),
                ["ATSC4i.json"] if descriptor["status"] == "complete" else [],
            )
        )
        result["descriptor"] = descriptor
        result["status"] = "complete"
        rows.sort(key=lambda r: ANALYSES.index(r["analysis"]))
        notes = [
            "Missing outputs remain explicit; equilibrium snapshots cannot supply surface pressure or tensile curves.",
            "Film geometry and zero guard-wall force must be verified alongside pressure statistics.",
            "Execution completion does not establish equilibration or convergence.",
        ]
        write_json(out / "analysis.json", result)
        report_markdown(out, "Part 3 analysis", rows, notes)
        print(out / "report.md")
    except Exception as e:
        result.update(status="failed", error=str(e))
        write_json(out / "analysis.json", result)
        raise


if __name__ == "__main__":
    main()
