#!/usr/bin/env python3
"""General spatial, segregation and component-resolved dynamics controller."""

import argparse
import collections
import json
import math
from pathlib import Path
import sys
from siliconelab_analysis_common import *

ANALYSES = [
    "component_density_profiles",
    "surface_enrichment",
    "network_profiles",
    "surface_core_comparisons",
    "phase_separation",
    "concentration_structure_factor",
    "oil_chain_aggregation",
    "MSD",
    "diffusion",
    "layer_dynamics",
    "Debye_Waller",
    "fixed_lag_dynamics",
]


def surfaces(profile, geometry, np):
    z = np.array(sorted(set(float(r["z_A"]) for r in profile)))
    total = np.array(
        [
            sum(float(r["density_g_cm3"]) for r in profile if float(r["z_A"]) == v)
            for v in z
        ]
    )
    smooth = np.convolve(total, np.ones(5) / 5, mode="same")
    peak = float(np.max(smooth))
    if geometry == "bulk":
        return {
            "status": "not_applicable",
            "reason": "periodic bulk has no material surfaces",
        }
    if peak <= 0:
        raise ValueError("Empty density profile")
    inside = np.flatnonzero(smooth >= 0.5 * peak)
    if len(inside) < 5:
        return {
            "status": "unresolved",
            "reason": "Too few material bins for half-maximum surfaces",
        }
    lo, hi = float(z[inside[0]]), float(z[inside[-1]])
    if inside[0] < 2 or inside[-1] >= len(z) - 2:
        return {
            "status": "unresolved",
            "reason": "Material touches the cell edge; no resolved vacuum interval",
        }
    if np.any(np.diff(inside) > 1):
        return {
            "status": "unresolved",
            "reason": "Disconnected half-maximum intervals; inspect film geometry",
        }
    return {
        "status": "complete",
        "z_lower_A": lo,
        "z_upper_A": hi,
        "thickness_A": hi - lo,
        "definition": "5-bin smoothed total mass density >= half its maximum; bin-center resolution",
    }


def enrichment(profile, surface, out):
    if surface["status"] != "complete":
        return None
    lo, hi = surface["z_lower_A"], surface["z_upper_A"]
    h = hi - lo
    regions = {
        "lower_surface": [lo, lo + 0.2 * h],
        "core": [lo + 0.2 * h, hi - 0.2 * h],
        "upper_surface": [hi - 0.2 * h, hi],
    }
    total_mass = sum(float(r["mass_g_per_mol"]) for r in profile)
    oil_mass = sum(
        float(r["mass_g_per_mol"]) for r in profile if r["component"] == "oil"
    )
    global_w = oil_mass / total_mass
    rows = []
    for name, (a, b) in regions.items():
        chosen = [r for r in profile if a <= float(r["z_A"]) < b]
        m = sum(float(r["mass_g_per_mol"]) for r in chosen)
        oil = sum(float(r["mass_g_per_mol"]) for r in chosen if r["component"] == "oil")
        rows.append(
            {
                "region": name,
                "z_lower_A": a,
                "z_upper_A": b,
                "mass_g_per_mol": m,
                "oil_mass_fraction": oil / m if m else None,
                "oil_enrichment_ratio": (
                    oil / m / global_w if m and global_w else None
                ),
            }
        )
    write_tsv(out / "surface_core.tsv", rows)
    return rows


def segregation(atoms, box, out, np):
    xyz = np.array([[float(a[k]) for k in ("x", "y", "z")] for a in atoms])
    lo = np.array(box["lo"])
    hi = np.array(box["hi"])
    lengths = hi - lo
    xyz[:, :2] = lo[:2] + np.mod(xyz[:, :2] - lo[:2], lengths[:2])
    if box["geometry"] == "bulk":
        xyz[:, 2] = lo[2] + np.mod(xyz[:, 2] - lo[2], lengths[2])
    sizes = np.maximum(1, np.ceil(lengths / 8).astype(int))
    indices = np.clip(((xyz - lo) / lengths * sizes).astype(int), 0, sizes - 1)
    flat = np.ravel_multi_index(indices.T, tuple(sizes))
    weight = np.array([float(a["mass"]) for a in atoms])
    selections = [
        (
            "oil_component",
            np.array([a["component"] == "oil" for a in atoms]),
            np.ones(len(atoms), dtype=bool),
            weight,
        ),
        (
            "MPS_chemistry",
            np.array([int(a["type"]) == 4 for a in atoms]),
            np.array([int(a["type"]) in (1, 2, 3, 4) for a in atoms]),
            np.ones(len(atoms)),
        ),
    ]
    result = {}
    for name, positive, valid, mass in selections:
        if not np.any(positive[valid]) or np.all(positive[valid]):
            result[name] = {
                "status": "not_applicable",
                "reason": "Only one population is present",
            }
            continue
        total = np.bincount(
            flat[valid], weights=mass[valid], minlength=int(np.prod(sizes))
        ).reshape(tuple(sizes))
        selected = np.bincount(
            flat[valid],
            weights=mass[valid] * positive[valid],
            minlength=int(np.prod(sizes)),
        ).reshape(tuple(sizes))
        w2 = np.bincount(
            flat[valid], weights=mass[valid] ** 2, minlength=int(np.prod(sizes))
        ).reshape(tuple(sizes))
        fraction = float(selected.sum() / total.sum())
        fields = {}
        for dim in ["3D", "XY"]:
            t = total if dim == "3D" else total.sum(axis=2)
            s = selected if dim == "3D" else selected.sum(axis=2)
            q = w2 if dim == "3D" else w2.sum(axis=2)
            occupied = t > 0
            f = np.divide(s, t, out=np.zeros_like(s), where=occupied)
            variance = float(np.sum(t * (f - fraction) ** 2) / t.sum())
            noise = float(
                fraction
                * (1 - fraction)
                * np.sum(np.divide(q, t, out=np.zeros_like(q), where=occupied))
                / t.sum()
            )
            excess = (variance - noise) / (fraction * (1 - fraction))
            # Raw FFT intensity of integrated cell concentration contrast; finite grid form factor not corrected.
            delta = s - fraction * t
            fft = np.fft.fftn(delta)
            intensity = np.abs(fft) ** 2 / float(np.sum(mass[valid] ** 2))
            axes = [
                2 * np.pi * np.fft.fftfreq(n, d=l / n)
                for n, l in zip(delta.shape, lengths[: delta.ndim])
            ]
            qgrid = np.sqrt(sum(g**2 for g in np.meshgrid(*axes, indexing="ij")))
            shell_width = float(2 * np.pi / max(lengths[: delta.ndim]))
            shell = np.rint(qgrid / shell_width).astype(int)
            rows = []
            maxq = min(
                0.8,
                float(np.min(np.pi / (lengths[: delta.ndim] / sizes[: delta.ndim]))),
            )
            for k in range(1, int(maxq / shell_width) + 1):
                mask = (shell == k) & (qgrid <= maxq)
                if mask.any():
                    rows.append(
                        {
                            "q_Ainv": float(qgrid[mask].mean()),
                            "Scc_grid_raw": float(intensity[mask].mean()),
                            "modes": int(mask.sum()),
                        }
                    )
            write_tsv(out / f"structure_factor_{name}_{dim}.tsv", rows)
            peak = max(rows, key=lambda r: r["Scc_grid_raw"]) if rows else None
            fields[dim] = {
                "weighted_composition_variance": variance,
                "independent_label_noise_estimate": noise,
                "excess_segregation_index": float(excess),
                "grid_dimensions": list(map(int, delta.shape)),
                "q_peak_Ainv": peak["q_Ainv"] if peak else None,
                "spacing_2pi_over_q_A": 2 * math.pi / peak["q_Ainv"] if peak else None,
                "peak_at_lowest_shell": bool(peak and peak == rows[0]),
                "spectrum": "raw coarse-grid concentration contrast; no particle form-factor correction",
            }
        result[name] = {
            "status": "complete",
            "global_fraction": fraction,
            "fraction_basis": (
                "mass" if name == "oil_component" else "backbone repeat count"
            ),
            "fields": fields,
        }
    write_json(out / "segregation.json", result)
    return result


def aggregation(atoms, box, out, np, cutoff):
    oil = [a for a in atoms if a["component"] == "oil"]
    molecules = sorted(set(int(a["molecule"]) for a in oil))
    if not molecules:
        return None
    parent = {m: m for m in molecules}
    contacts = set()
    cells = collections.defaultdict(list)
    lo = np.array(box["lo"])
    lengths = np.array(box["hi"]) - lo
    n = np.maximum(1, np.floor(lengths / cutoff).astype(int))
    xyz = []
    for a in oil:
        p = np.array([float(a[k]) for k in ("x", "y", "z")])
        p[:2] = lo[:2] + np.mod(p[:2] - lo[:2], lengths[:2])
        if box["geometry"] == "bulk":
            p[2] = lo[2] + np.mod(p[2] - lo[2], lengths[2])
        idx = tuple(np.clip(((p - lo) / lengths * n).astype(int), 0, n - 1))
        xyz.append(p)
        cells[idx].append(len(xyz) - 1)

    def find(m):
        while parent[m] != m:
            parent[m] = parent[parent[m]]
            m = parent[m]
        return m

    for cell, ids in cells.items():
        neighbors = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    q = [cell[0] + dx, cell[1] + dy, cell[2] + dz]
                    q[0] %= int(n[0])
                    q[1] %= int(n[1])
                    if box["geometry"] == "bulk":
                        q[2] %= int(n[2])
                    if 0 <= q[2] < n[2]:
                        neighbors.add(tuple(q))
        for a in ids:
            for q in neighbors:
                for b in cells.get(q, []):
                    if b <= a:
                        continue
                    ma, mb = int(oil[a]["molecule"]), int(oil[b]["molecule"])
                    if ma == mb:
                        continue
                    diff = xyz[a] - xyz[b]
                    diff[:2] -= np.rint(diff[:2] / lengths[:2]) * lengths[:2]
                    if box["geometry"] == "bulk":
                        diff[2] -= round(diff[2] / lengths[2]) * lengths[2]
                    if float(diff @ diff) < cutoff**2:
                        contacts.add(tuple(sorted((ma, mb))))
                        ra, rb = find(ma), find(mb)
                        parent[ra] = rb
    groups = collections.Counter(find(m) for m in molecules)
    r = {
        "status": "complete",
        "definition": "Generic oil-oil bead contact clusters; includes pendant beads; snapshot only",
        "contact_cutoff_A": cutoff,
        "oil_chains": len(molecules),
        "clusters": len(groups),
        "largest_cluster_fraction": max(groups.values()) / len(molecules),
        "contacting_molecule_pairs": len(contacts),
        "warning": "Operational contact connectivity is not proof of thermodynamic aggregation; vary cutoff.",
    }
    write_json(out / "oil_clusters.json", r)
    return r


def summarize_dynamics(folder, out, kind, np):
    sampling = json.loads((folder / "sampling.json").read_text())
    rows = read_tsv(folder / "dynamics.tsv")
    globalrows = [r for r in rows if r["layer_center_relative_A"] == "global"]
    duration = (
        (sampling["last_step"] - sampling["first_step"])
        * sampling["timestep_fs"]
        * 0.001
    )
    fixed = 10 if kind == "dw" else 10000
    fixedrows = [dict(r) for r in rows if abs(float(r["lag_ps"]) - fixed) < 1e-8]
    if kind == "dw":
        globalvalues = {
            (r["component"], r["particle"]): r
            for r in fixedrows
            if r["layer_center_relative_A"] == "global"
        }
        for r in fixedrows:
            g = globalvalues[(r["component"], r["particle"])]
            for component in ["xy", "3d"]:
                u = float(r["MSD_" + component + "_A2"])
                base = float(g["MSD_" + component + "_A2"])
                r["inverse_u2_" + component + "_Ainv2"] = 1 / u if u > 0 else None
                r["u2_" + component + "_over_global"] = u / base if base > 0 else None
    write_tsv(out / ("fixed_lag_" + kind + ".tsv"), fixedrows)
    fits = []
    for key in sorted(set((r["component"], r["particle"]) for r in globalrows)):
        selected = [
            r
            for r in globalrows
            if (r["component"], r["particle"]) == key
            and 0.5 * duration <= float(r["lag_ps"]) <= 0.9 * duration
        ]
        if len(selected) < 3:
            continue
        x = np.array([float(r["lag_ps"]) for r in selected])
        y = np.array([float(r["MSD_xy_A2"]) for r in selected])
        coef = np.polyfit(x, y, 1)
        res = y - np.polyval(coef, x)
        sst = float(np.sum((y - y.mean()) ** 2))
        r2 = 1 - float(np.sum(res**2)) / sst if sst else None
        directional = {}
        for dimension, factor in [("x", 2), ("y", 2), ("z", 2), ("xy", 4), ("3d", 6)]:
            values = np.array([float(r["MSD_" + dimension + "_A2"]) for r in selected])
            slope = float(np.polyfit(x, values, 1)[0])
            directional["D_" + dimension + "_A2_ps"] = slope / factor
        fits.append(
            dict(
                directional,
                component=key[0],
                particle=key[1],
                fit_start_ps=float(x.min()),
                fit_end_ps=float(x.max()),
                R2_xy=r2,
                status="preliminary_fit" if coef[0] >= 0 else "invalid_negative_slope",
                interpretation=(
                    "Network motion/relaxation"
                    if key[0] != "oil"
                    else (
                        "Oil translation"
                        if key[1] == "molecular_COM"
                        else "Oil bead motion includes internal relaxation"
                    )
                ),
                film_warning="Confined z and 3D slopes are comparison diagnostics, not asymptotic self-diffusion",
            )
        )
    if kind == "msd":
        write_json(out / "diffusion.json", fits)
    return {
        "sampling": sampling,
        "duration_ps": duration,
        "fixed_lag_ps": fixed,
        "fixed_lag_available": bool(fixedrows),
        "fits": fits,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data", type=Path)
    p.add_argument("info", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--geometry", choices=["bulk", "film"])
    p.add_argument("--trajectory", type=Path)
    p.add_argument("--dw-trajectory", type=Path)
    p.add_argument("--backend-dir")
    p.add_argument("--origin-stride", type=int, default=100)
    p.add_argument("--max-origins", type=int, default=16)
    p.add_argument("--contact-cutoff", type=float, default=8)
    a = p.parse_args(argv)
    try:
        import numpy as np
    except ImportError:
        raise RuntimeError("Part2 requires NumPy; see requirements-analysis.txt")
    if a.output_dir.exists():
        raise ValueError("Choose a new output directory; results are never overwritten")
    if a.contact_cutoff <= 0 or a.origin_stride < 1 or a.max_origins < 1:
        raise ValueError("Positive sampling parameters required")
    source, normalized = metadata(a.info, a.geometry)
    out = a.output_dir
    out.mkdir(parents=True)
    write_json(out / "normalized.info", normalized)
    result = {
        "status": "running",
        "source": {
            "data": str(a.data.resolve()),
            "data_sha256": sha(a.data),
            "info": str(a.info.resolve()),
            "info_sha256": sha(a.info),
        },
        "analyses": [],
    }
    write_json(out / "analysis.json", result)
    try:
        invoke(
            [
                backend("siliconelab_static", a.backend_dir),
                str(a.data.resolve()),
                str((out / "normalized.info").resolve()),
                str((out / "validation").resolve()),
            ],
            out,
        )
        invoke(
            [
                backend("siliconelab_spatial_dynamics", a.backend_dir),
                "spatial",
                str(a.data.resolve()),
                str((out / "normalized.info").resolve()),
                str(out.resolve()),
            ],
            out,
        )
        profiles = read_tsv(out / "component_profiles.tsv")
        atoms = read_tsv(out / "snapshot_atoms.tsv")
        box = json.loads((out / "box.json").read_text())
        surface = surfaces(profiles, normalized["geometry"], np)
        write_json(out / "surfaces.json", surface)
        enrichment(profiles, surface, out)
        mixed = bool(
            normalized["components"]["filler"]["M"] and normalized["network_enabled"]
        )
        rows = result["analyses"]
        rows.append(
            state_row(ANALYSES[0], "complete", files=["component_profiles.tsv"])
        )
        rows.append(
            state_row(
                ANALYSES[1],
                (
                    (
                        "complete"
                        if surface["status"] == "complete"
                        else "insufficient_geometry"
                    )
                    if mixed and normalized["geometry"] == "film"
                    else "not_applicable"
                ),
                "Oil/network mixture and resolved film surfaces required",
                (
                    ["surface_core.tsv"]
                    if mixed and surface["status"] == "complete"
                    else []
                ),
            )
        )
        if normalized["network_enabled"]:
            invoke(
                [
                    backend("siliconelab_profiles", a.backend_dir),
                    str(a.data.resolve()),
                    str((out / "normalized.info").resolve()),
                    "--output-dir",
                    str((out / "network").resolve()),
                    "--no-z1",
                ],
                out,
            )
            rows.append(
                state_row(
                    ANALYSES[2],
                    "complete",
                    "Native local conversion, defects and strand geometry; native wall-coordinate summaries are auxiliary for free films",
                    ["network/"],
                )
            )
        else:
            rows.append(state_row(ANALYSES[2], "not_applicable", "No reactive network"))
        rows.append(
            state_row(
                ANALYSES[3],
                (
                    "complete"
                    if surface["status"] == "complete"
                    else (
                        "not_applicable"
                        if normalized["geometry"] == "bulk"
                        else "insufficient_geometry"
                    )
                ),
                "Measured density surfaces, 20% edge regions and 60% core",
                (
                    ["surfaces.json", "surface_core.tsv"]
                    if surface["status"] == "complete"
                    else []
                ),
            )
        )
        if (
            source.get("format") in ("V22-model-info", "V35-model-info")
            and source.get("silicone_oil", {}).get("total_mps_repeats", 0) > 0
        ):
            phase = out / "chemical_phase_native"
            phase.mkdir()
            invoke(
                [
                    backend("siliconelab_phase", a.backend_dir),
                    str(a.data.resolve()),
                    str(a.info.resolve()),
                    "--output",
                    str((phase / "metrics.dat").resolve()),
                    "--structure-output",
                    str((phase / "structure_factor.dat").resolve()),
                    "--field-output",
                    str((phase / "composition_field.dat").resolve()),
                    "--report-output",
                    str((phase / "report.txt").resolve()),
                ],
                out,
            )
        seg = segregation(atoms, box, out, np)
        available = any(v["status"] == "complete" for v in seg.values())
        rows.extend(
            [
                state_row(
                    ANALYSES[i],
                    "complete" if available else "not_applicable",
                    "Chemical and component fields are reported separately; raw grid spectra require grid/cutoff/box checks",
                    ["segregation.json"],
                )
                for i in (4, 5)
            ]
        )
        cluster = aggregation(atoms, box, out, np, a.contact_cutoff)
        rows.append(
            state_row(
                ANALYSES[6],
                "complete" if cluster else "not_applicable",
                files=["oil_clusters.json"] if cluster else [],
            )
        )
        dyn = {}
        for kind, path in [("msd", a.trajectory), ("dw", a.dw_trajectory)]:
            if not path:
                continue
            if not path.is_file():
                dyn[kind] = {"status": "missing_output"}
                continue
            folder = out / kind
            folder.mkdir()
            invoke(
                [
                    backend("siliconelab_spatial_dynamics", a.backend_dir),
                    "dynamics",
                    str(a.data.resolve()),
                    str((out / "normalized.info").resolve()),
                    str(folder.resolve()),
                    str(path.resolve()),
                    str(a.origin_stride),
                    str(a.max_origins),
                ],
                out,
            )
            dyn[kind] = summarize_dynamics(folder, out, kind, np)
            dyn[kind]["status"] = "complete"
            dyn[kind]["source_sha256"] = sha(path)
        has_msd = dyn.get("msd", {}).get("status") == "complete"
        has_dw = dyn.get("dw", {}).get("status") == "complete"
        rows.append(
            state_row(
                ANALYSES[7],
                "complete" if has_msd else "missing_output",
                files=["msd/dynamics.tsv"] if has_msd else [],
            )
        )
        rows.append(
            state_row(
                ANALYSES[8],
                "preliminary" if has_msd else "missing_output",
                "Fit alone does not establish long-time diffusion; inspect oil COM and fit interval",
                ["diffusion.json"] if has_msd else [],
            )
        )
        rows.append(
            state_row(
                ANALYSES[9],
                "complete" if has_msd else "missing_output",
                "Component bead layers use each origin and material COM reference",
                ["msd/dynamics.tsv"] if has_msd else [],
            )
        )
        dwgood = has_dw and dyn["dw"]["fixed_lag_available"]
        rows.append(
            state_row(
                ANALYSES[10],
                (
                    "complete"
                    if dwgood
                    else ("insufficient_sampling" if has_dw else "missing_output")
                ),
                "Explicit 10 ps displacement; no automatic plateau claim; 1/u2 is a proxy",
                ["dw/dynamics.tsv", "fixed_lag_dw.tsv"] if dwgood else [],
            )
        )
        rows.append(
            state_row(
                ANALYSES[11],
                (
                    "complete"
                    if has_msd and dyn["msd"]["fixed_lag_available"] and dwgood
                    else "insufficient_sampling"
                ),
                "Requires exact 10 ns and 10 ps partner frames",
                ["fixed_lag_msd.tsv", "fixed_lag_dw.tsv"],
            )
        )
        result.update(
            status="complete",
            geometry=normalized["geometry"],
            surface=surface,
            dynamics=dyn,
            numpy_version=np.__version__,
            notes=[
                "Static profiles and segregation are snapshot measurements; trajectory averaging/coarsening requires multiple analyzed snapshots.",
                "Sampled time origins overlap and are correlated; no independent-replica error is inferred.",
                "Density-derived film boundaries are operational; inspect profiles and vacuum.",
                "Generic oil contacts and MPS-specific aggregation must not be conflated.",
                "Z1 spatial profiles need a validated contour mapping; mixed network/oil exports are not blindly imported.",
            ],
        )
        if (
            sha(a.data) != result["source"]["data_sha256"]
            or sha(a.info) != result["source"]["info_sha256"]
        ):
            raise RuntimeError("Inputs changed during analysis")
        write_json(out / "analysis.json", result)
        report_markdown(out, "Part 2 analysis", rows, result["notes"])
        print(out / "report.md")
    except Exception as e:
        result.update(status="failed", error=str(e))
        write_json(out / "analysis.json", result)
        raise


if __name__ == "__main__":
    main()
