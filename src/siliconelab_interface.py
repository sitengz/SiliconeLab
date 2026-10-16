"""Bulk rigid-NP interface analysis with explicit configuration and applicability.

Distances are from the nominal sphere of constituent bead centers, reconstructed
from each current snapshot. Native chemical graph labels are not NP contacts or
primitive-path entanglements. All calculations preserve the original inputs.
"""

import argparse
import collections
import json
import math
from pathlib import Path
import sys

import numpy as np
from siliconelab_analysis_common import (
    backend,
    invoke,
    read_tsv,
    sha,
    static_driver,
    write_json,
    write_tsv,
)

DEFAULTS = dict(
    shell_edges_A=[0, 5, 10, 15, 20, 30, 40, 60],
    interface_gap_A=10.0,
    reference_gap_A=30.0,
    volume_points=250000,
    volume_replicates=4,
    seed=20261016,
    max_dynamic_atoms=20000,
    max_origins=32,
    origin_stride=1,
    dw_lags_ps=[0.1, 0.5, 1.0, 4.0, 10.0, 25.0, 50.0],
    msd_lags_ps=[25.0, 100.0, 500.0, 1000.0, 2500.0, 5000.0, 10000.0],
    timestep_fs=5.0,
    max_frame_displacement_A=20.0,
    rigidity_tolerance_A=0.02,
    tensile_fit_windows=[[0.005, 0.02], [0.01, 0.03], [0.02, 0.05]],
    tensile_strain_block=0.002,
    tensile_profile_interval_strain=0.1,
    comparison_windows_A=[[0, 5], [5, 10], [0, 10], [5, 15], [10, 20]],
    comparison_reference_gaps_A=[20, 30, 40],
    plot=True,
)
PATHS = {"data", "info", "output_dir", "backend_dir", "part2_report"}
EXTRA = {"trajectories", "tensile_tables", "tensile_data", "tensile_trajectories"}


def configure(path):
    p = Path(path).resolve()
    raw = json.loads(p.read_text())
    if set(raw) - set(DEFAULTS) - PATHS - EXTRA:
        raise ValueError(
            "Unknown interface configuration keys: "
            + str(set(raw) - set(DEFAULTS) - PATHS - EXTRA)
        )
    if not {"data", "info", "output_dir"}.issubset(raw):
        raise ValueError("data, info and output_dir required")
    c = {**DEFAULTS, **raw}
    for k in PATHS:
        if c.get(k):
            c[k] = str((p.parent / str(c[k])).resolve())
    c["trajectories"] = [
        {**t, "path": str((p.parent / t["path"]).resolve())}
        for t in c.get("trajectories", [])
    ]
    if any(
        set(t) - {"path", "kind", "label", "timestep_fs"}
        or t.get("kind") not in ("dw", "msd")
        for t in c["trajectories"]
    ):
        raise ValueError("Trajectory entries require path and kind dw/msd")
    for key in ("tensile_tables", "tensile_data", "tensile_trajectories"):
        c[key] = {a: str((p.parent / s).resolve()) for a, s in c.get(key, {}).items()}
        if set(c[key]) - set("xyz"):
            raise ValueError("Tensile axes must be x/y/z")
    for key in (
        "interface_gap_A",
        "reference_gap_A",
        "timestep_fs",
        "max_frame_displacement_A",
        "rigidity_tolerance_A",
        "tensile_strain_block",
        "tensile_profile_interval_strain",
    ):
        if (
            isinstance(c[key], bool)
            or not isinstance(c[key], (int, float))
            or not math.isfinite(c[key])
            or c[key] <= 0
        ):
            raise ValueError("Positive finite " + key + " required")
    for key in (
        "volume_points",
        "volume_replicates",
        "max_dynamic_atoms",
        "max_origins",
        "origin_stride",
        "seed",
    ):
        if type(c[key]) is not int or c[key] <= 0:
            raise ValueError("Positive integer " + key + " required")
    edges = c["shell_edges_A"]
    if (
        not edges
        or edges[0] != 0
        or any(not math.isfinite(v) or b <= v for v, b in zip(edges, edges[1:]))
        or not math.isfinite(edges[-1])
    ):
        raise ValueError("Shell edges must start at zero and strictly increase")
    if c["reference_gap_A"] <= c["interface_gap_A"] or any(
        c[k] not in edges for k in ("interface_gap_A", "reference_gap_A")
    ):
        raise ValueError("Interface/reference boundaries must be separate shell edges")
    for key in ("dw_lags_ps", "msd_lags_ps"):
        if not c[key] or any(not math.isfinite(v) or v <= 0 for v in c[key]):
            raise ValueError("Positive finite lag times required")
    for lo, hi in c["tensile_fit_windows"]:
        if not 0 <= lo < hi or not math.isfinite(hi):
            raise ValueError("Invalid tensile fit window")
    for lo, hi in c["comparison_windows_A"]:
        if not 0 <= lo < hi or not math.isfinite(hi):
            raise ValueError("Invalid interface comparison window")
    if any(not math.isfinite(x) or x <= 0 for x in c["comparison_reference_gaps_A"]):
        raise ValueError("Invalid comparison reference gap")
    if type(c["plot"]) is not bool:
        raise ValueError("plot must be a boolean")
    c["config_path"] = str(p)
    return c


def mic(x, L):
    return x - L * np.rint(x / L)


def snapshot(path, expected):
    pos = np.empty((expected, 3))
    mol = np.zeros(expected, dtype=np.int64)
    typ = np.zeros(expected, dtype=np.int32)
    bonds = []
    bounds = {}
    masses = {}
    section = None
    count = 0
    seen = set()
    declared = None
    with Path(path).open() as f:
        next(f)
        for line in f:
            w = line.split("#", 1)[0].split()
            if not w:
                continue
            if len(w) == 2 and w[1] == "atoms":
                declared = int(w[0])
            if len(w) == 4 and w[2] in ("xlo", "ylo", "zlo"):
                bounds[w[2][0]] = tuple(map(float, w[:2]))
            if not w[0].lstrip("-").isdigit():
                section = " ".join(w)
                continue
            if section == "Atoms":
                i, m, t = map(int, w[:3])
                i -= 1
                if i in seen or not 0 <= i < expected or len(w) not in (7, 10):
                    raise ValueError("Invalid atom identity/record")
                seen.add(i)
                mol[i] = m
                typ[i] = t
                pos[i] = list(map(float, w[4:7]))
                count += 1
            elif section == "Bonds":
                bonds.append(tuple(map(int, w)))
            elif section == "Masses":
                masses[int(w[0])] = float(w[1])
    if (
        declared != expected
        or count != expected
        or len(bounds) != 3
        or not np.isfinite(pos).all()
    ):
        raise ValueError("Snapshot count/coordinate contract failed")
    low = np.array([bounds[a][0] for a in "xyz"])
    L = np.array([bounds[a][1] - bounds[a][0] for a in "xyz"])
    if not np.isfinite(L).all() or (L <= 0).any():
        raise ValueError("Invalid box")
    return dict(
        pos=pos,
        mol=mol,
        typ=typ,
        bonds=bonds,
        L=L,
        low=low,
        mass=np.array([masses[int(t)] for t in typ]),
    )


def body_centers(pos, mol, info, L):
    centers = []
    reference = []
    if 4 * info["np_radius_bead_centers_A"] >= min(L):
        raise ValueError(
            "NP diameter must be below half every box length for local body reconstruction"
        )
    for m in range(info["np_molecule_id_start"], info["total_molecules"] + 1):
        ids = np.flatnonzero(mol == m)
        if len(ids) != info["beads_per_np"]:
            raise ValueError("Incomplete NP body")
        xyz = pos[ids[0]] + mic(pos[ids] - pos[ids[0]], L)
        centers.append(xyz.mean(axis=0))
        reference.append(np.linalg.norm(xyz[:4, None, :] - xyz[None, :, :], axis=2))
    return np.array(centers), reference


def nearest(pos, centers, L, radius):
    gaps = np.empty(len(pos))
    which = np.empty(len(pos), dtype=int)
    normal = np.empty((len(pos), 3))
    for start in range(0, len(pos), 4096):
        delta = mic(pos[start : start + 4096, None, :] - centers[None, :, :], L)
        distance = np.linalg.norm(delta, axis=2)
        n = distance.argmin(axis=1)
        v = delta[np.arange(len(n)), n]
        r = distance[np.arange(len(n)), n]
        gaps[start : start + len(n)] = r - radius
        which[start : start + len(n)] = n
        normal[start : start + len(n)] = np.divide(
            v, r[:, None], out=np.zeros_like(v), where=r[:, None] > 0
        )
    return gaps, which, normal


def shell_index(gap, c):
    return np.searchsorted(c["shell_edges_A"], gap, side="right")


def shell_record(i, c):
    edges = c["shell_edges_A"]
    return dict(
        shell=i,
        gap_lower_A=edges[i - 1] if i else None,
        gap_upper_A=edges[i] if i < len(edges) else None,
        region=(
            "nominal_NP_interior"
            if i == 0
            else (
                (
                    "interface"
                    if edges[i] <= c["interface_gap_A"]
                    else (
                        "reference"
                        if edges[i - 1] >= c["reference_gap_A"]
                        else "transition"
                    )
                )
                if i < len(edges)
                else "reference"
            )
        ),
    )


def volumes(s, centers, info, c):
    rng = np.random.default_rng(c["seed"])
    nb = len(c["shell_edges_A"]) + 1
    counts = []
    per = []
    V = float(np.prod(s["L"]))
    for _ in range(c["volume_replicates"]):
        points = s["low"] + rng.random((c["volume_points"], 3)) * s["L"]
        gap, which, _ = nearest(
            points, centers, s["L"], info["np_radius_bead_centers_A"]
        )
        bins = shell_index(gap, c)
        counts.append(np.bincount(bins, minlength=nb) * V / len(points))
        per.append(
            np.bincount(which * nb + bins, minlength=len(centers) * nb).reshape(
                len(centers), nb
            )
            * V
            / len(points)
        )
    a = np.array(counts)
    return (
        a.mean(axis=0),
        (
            a.std(axis=0, ddof=1) / math.sqrt(len(a))
            if len(a) > 1
            else np.full(nb, np.nan)
        ),
        np.mean(per, axis=0),
    )


def pruned_core(contours):
    active = [s for s in contours if s["status"] == "active"]
    adj = collections.defaultdict(list)
    for s in active:
        for u, v in ((s["u"], s["v"]), (s["v"], s["u"])):
            adj[u].append((v, s["id"]))
    components = []
    seen = set()
    for root in adj:
        if root in seen:
            continue
        stack = [root]
        nodes = set()
        edges = set()
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            nodes.add(u)
            for v, e in adj[u]:
                edges.add(e)
                stack.append(v)
        components.append((nodes, edges))
    if not components:
        return set(), set()
    nodes, edges = max(components, key=lambda r: (len(r[1]), len(r[0])))
    degree = {u: sum(e in edges for _, e in adj[u]) for u in nodes}
    q = collections.deque(u for u in nodes if degree[u] < 2)
    core = set(edges)
    while q:
        u = q.popleft()
        for v, e in adj[u]:
            if e not in core:
                continue
            core.remove(e)
            degree[v] -= 1
            if degree[v] == 1:
                q.append(v)
    return core, edges


def network_map(folder, s, normalized):
    contours = []
    for r in read_tsv(folder / "contours.tsv"):
        contours.append(
            dict(
                id=int(r["strand_id"]),
                parent=int(r["parent_molecule"]),
                status=r["status"],
                u=int(r["first_node"]),
                v=int(r["second_node"]),
                ids=np.fromstring(r["atom_ids"], sep=",", dtype=np.int64) - 1,
            )
        )
    nodes = read_tsv(folder / "nodes.tsv")
    core, giant = pruned_core(contours)
    n = len(s["pos"])
    weights = {
        k: np.zeros(n)
        for k in (
            "all",
            "active",
            "core",
            "dangling",
            "self_loop",
            "dangling_loop",
            "isolated",
        )
    }
    masks = {
        k: np.zeros(n, dtype=bool)
        for k in (
            "native_active",
            "two_core_candidate",
            "defective_strand",
            "strand_midpoint",
            "junction_attachment",
        )
    }
    nodepos = {
        int(r["node_id"]): np.array([float(r[a + "_A"]) for a in "xyz"]) for r in nodes
    }
    adj = collections.defaultdict(list)
    for strand in contours:
        ids = strand["ids"]
        strand["two_core_candidate"] = strand["id"] in core
        if (
            (ids < 0).any()
            or (ids >= normalized["composition"]["total_beads"]).any()
            or not (s["mol"][ids] == strand["parent"]).all()
        ):
            raise ValueError("Native contour ownership mismatch")
        key = strand["status"]
        np.add.at(weights["all"], ids, 1 / len(ids))
        np.add.at(weights[key], ids, 1 / len(ids))
        if strand["id"] in core:
            np.add.at(weights["core"], ids, 1 / len(ids))
            masks["two_core_candidate"][ids] = True
        masks["native_active" if key == "active" else "defective_strand"][ids] = True
        masks["strand_midpoint"][ids[len(ids) // 2]] = True
        if key == "active":
            u, v = strand["u"], strand["v"]
            a = s["pos"][ids]
            d = (
                mic(a[0] - nodepos[u], s["L"])
                + mic(np.diff(a, axis=0), s["L"]).sum(axis=0)
                + mic(nodepos[v] - a[-1], s["L"])
            )
            shift = np.rint((d - (nodepos[v] - nodepos[u])) / s["L"]).astype(int)
            if np.max(np.abs(d - (nodepos[v] - nodepos[u]) - shift * s["L"])) > 1e-6:
                raise ValueError("Edge lattice shift is not integral")
            adj[u].append((v, shift))
            adj[v].append((u, -shift))
    for r in nodes:
        ids = np.fromstring(r["attachment_atoms"], sep=",", dtype=int) - 1
        masks["junction_attachment"][ids] = True
    potentials = {}
    cycles = []
    for root in adj:
        if root in potentials:
            continue
        potentials[root] = np.zeros(3, dtype=int)
        stack = [root]
        while stack:
            u = stack.pop()
            for v, shift in adj[u]:
                if v not in potentials:
                    potentials[v] = potentials[u] + shift
                    stack.append(v)
                else:
                    winding = potentials[u] + shift - potentials[v]
                    if winding.any():
                        cycles.append(winding.tolist())
    summary = dict(
        native_strands=len(contours),
        native_status_counts=dict(collections.Counter(s["status"] for s in contours)),
        largest_active_component_edges=len(giant),
        giant_component_two_core_candidate_edges=len(core),
        corrected_periodic_cycle_winding_rank=(
            int(np.linalg.matrix_rank(np.array(cycles))) if cycles else 0
        ),
        definitions={
            "native_active": "native effective strand connects two distinct covalent junctions",
            "two_core_candidate": "native-active edges remaining after degree<2 pruning in the largest active component; candidate backbone, not exact elastically active strand count",
            "periodic_winding": "bond-by-bond path plus junction attachments, relative to raw wrapped node separation; no inherited directional self-path convention",
        },
        caution="A 2-core can retain pendant cycles. These graph labels do not include NP-induced entanglements or prove mechanical activity.",
    )
    return contours, nodes, weights, masks, summary


def ratio(a, b):
    return (
        float(a / b)
        if a is not None
        and b is not None
        and b != 0
        and math.isfinite(a)
        and math.isfinite(b)
        else None
    )


def profile_snapshot(s, info, normalized, contours, nodes, weights, c, out):
    out.mkdir()
    centers, _ = body_centers(s["pos"], s["mol"], info, s["L"])
    matrix = s["typ"] != info["np_atom_type"]
    gap, particle, normal = nearest(
        s["pos"], centers, s["L"], info["np_radius_bead_centers_A"]
    )
    bins = shell_index(gap, c)
    volume, volume_se, per_volume = volumes(s, centers, info, c)
    nb = len(volume)
    histogram = lambda values: np.bincount(bins, weights=values, minlength=nb)
    count = histogram(matrix)
    mass = histogram(s["mass"] * matrix)
    component_counts = {}
    for name, comp in normalized["components"].items():
        chosen = (s["mol"] >= (comp["molecule_id_start"] or -1)) & (
            s["mol"] <= (comp["molecule_id_end"] or -1)
        )
        component_counts[name] = histogram(chosen)
    eq = {k: histogram(w) for k, w in weights.items()}
    node_gap, _, _ = nearest(
        np.array([[float(r[a + "_A"]) for a in "xyz"] for r in nodes]).reshape(-1, 3),
        centers,
        s["L"],
        info["np_radius_bead_centers_A"],
    )
    node_bins = shell_index(node_gap, c)
    node_count = np.bincount(node_bins, minlength=nb)
    degree = np.array([int(r["active_degree"]) for r in nodes])
    degree_sum = np.bincount(node_bins, weights=degree, minlength=nb)
    poor = np.bincount(node_bins, weights=(degree < 3), minlength=nb)
    reactive = {}
    reacted = set()
    bond_type = normalized["force_field"]["bond_type_map"]["crosslink"]
    for _, t, i, j in s["bonds"]:
        if t == bond_type:
            reacted.update((i - 1, j - 1))
    for name, sites in [
        ("strands", normalized["strand"]["reactive_bead_sites"]),
        ("crosslinkers", normalized["crosslinker"]["reactive_bead_sites"]),
    ]:
        comp = normalized["components"][name]
        ids = []
        for mol in range(
            comp["molecule_id_start"] or 0, (comp["molecule_id_end"] or -1) + 1
        ):
            members = np.flatnonzero(s["mol"] == mol)
            ids.extend(members[np.array(sites, dtype=int) - 1])
        ids = np.array(ids, dtype=int)
        possible = np.bincount(bins[ids], minlength=nb)
        taken = np.bincount(bins[[i for i in ids if i in reacted]], minlength=nb)
        reactive[name] = (possible, taken)
    orientation = np.zeros(nb)
    orientation_n = np.zeros(nb)
    strand_rows = []
    for strand in contours:
        ids = strand["ids"]
        a = s["pos"][ids]
        delta = mic(a[:, None, :] - centers[None, :, :], s["L"])
        contact = (
            np.linalg.norm(delta, axis=2) - info["np_radius_bead_centers_A"]
            <= c["interface_gap_A"]
        ).any(axis=0)
        path = np.vstack(
            [a[0], a[0] + np.cumsum(mic(np.diff(a, axis=0), s["L"]), axis=0)]
        )
        strand_rows.append(
            dict(
                strand_id=strand["id"],
                parent_molecule=strand["parent"],
                status=strand["status"],
                two_core_candidate=strand["two_core_candidate"],
                minimum_surface_gap_A=float(gap[ids].min()),
                mean_surface_gap_A=float(gap[ids].mean()),
                contacted_NPs=int(contact.sum()),
                proximity_class=(
                    "bridges_multiple_NPs"
                    if contact.sum() > 1
                    else "single_NP_proximity" if contact.any() else "no_NP_proximity"
                ),
                contour_length_A=float(
                    np.linalg.norm(np.diff(path, axis=0), axis=1).sum()
                ),
                Ree_A=float(np.linalg.norm(path[-1] - path[0])),
                Rg_A=float(
                    np.sqrt(np.mean(np.sum((path - path.mean(axis=0)) ** 2, axis=1)))
                ),
            )
        )
        if strand["status"] == "active":
            tangent = mic(np.diff(a, axis=0), s["L"])
            midpoint = a[:-1] + 0.5 * tangent
            g, _, n = nearest(
                midpoint, centers, s["L"], info["np_radius_bead_centers_A"]
            )
            b = shell_index(g, c)
            length = np.linalg.norm(tangent, axis=1)
            cos = np.divide(
                np.einsum("ij,ij->i", tangent, n),
                length,
                out=np.zeros_like(length),
                where=length > 0,
            )
            orientation += np.bincount(
                b, weights=0.5 * (3 * cos * cos - 1), minlength=nb
            )
            orientation_n += np.bincount(b, minlength=nb)
    rows = []
    for i, V in enumerate(volume):
        row = {
            **shell_record(i, c),
            "volume_A3": float(V),
            "volume_MC_SE_A3": (
                float(volume_se[i]) if math.isfinite(volume_se[i]) else None
            ),
            "matrix_beads": int(count[i]),
            "matrix_number_density_Ainv3": ratio(count[i], V),
            "matrix_mass_density_g_cm3": ratio(mass[i], 0.602214076 * V),
            "all_strand_equivalents": float(eq["all"][i]),
            "native_active_strand_equivalents": float(eq["active"][i]),
            "two_core_candidate_strand_equivalents": float(eq["core"][i]),
            "native_active_strand_equivalent_density_Ainv3": ratio(eq["active"][i], V),
            "two_core_candidate_strand_equivalent_density_Ainv3": ratio(
                eq["core"][i], V
            ),
            "native_active_fraction_of_strand_equivalents": ratio(
                eq["active"][i], eq["all"][i]
            ),
            "two_core_candidate_fraction_of_strand_equivalents": ratio(
                eq["core"][i], eq["all"][i]
            ),
            "dangling_fraction_of_strand_equivalents": ratio(
                eq["dangling"][i], eq["all"][i]
            ),
            "self_loop_fraction_of_strand_equivalents": ratio(
                eq["self_loop"][i], eq["all"][i]
            ),
            "isolated_fraction_of_strand_equivalents": ratio(
                eq["isolated"][i], eq["all"][i]
            ),
            "junctions": int(node_count[i]),
            "junction_density_Ainv3": ratio(node_count[i], V),
            "mean_active_junction_degree": ratio(degree_sum[i], node_count[i]),
            "junction_degree_below_three_fraction": ratio(poor[i], node_count[i]),
            "native_active_bond_radial_P2": ratio(orientation[i], orientation_n[i]),
        }
        for name, (possible, taken) in reactive.items():
            row[name + "_site_conversion"] = ratio(taken[i], possible[i])
            row[name + "_possible_sites"] = int(possible[i])
            row[name + "_reacted_sites"] = int(taken[i])
        for name, values in component_counts.items():
            row[name + "_beads"] = int(values[i])
            row[name + "_fraction_of_matrix_beads"] = ratio(values[i], count[i])
        rows.append(row)
    write_tsv(out / "radial_profiles.tsv", rows)
    write_tsv(out / "strand_proximity.tsv", strand_rows)
    per = []
    for p in range(len(centers)):
        chosen = particle == p
        nc = np.bincount(bins, weights=matrix * chosen, minlength=nb)
        ac = np.bincount(bins, weights=weights["active"] * chosen, minlength=nb)
        allc = np.bincount(bins, weights=weights["all"] * chosen, minlength=nb)
        for i in range(nb):
            per.append(
                dict(
                    np_molecule=p + info["np_molecule_id_start"],
                    **shell_record(i, c),
                    volume_A3=float(per_volume[p, i]),
                    matrix_beads=int(nc[i]),
                    matrix_number_density_Ainv3=ratio(nc[i], per_volume[p, i]),
                    native_active_strand_equivalents=float(ac[i]),
                    all_strand_equivalents=float(allc[i]),
                    native_active_fraction=ratio(ac[i], allc[i]),
                )
            )
    write_tsv(out / "per_particle_profiles.tsv", per)
    summary = {}
    for region in ("interface", "reference"):
        selected = [i for i, r in enumerate(rows) if r["region"] == region]
        V = sum(volume[i] for i in selected)
        summary[region] = dict(
            volume_A3=float(V),
            matrix_beads=int(sum(count[i] for i in selected)),
            matrix_number_density_Ainv3=ratio(sum(count[i] for i in selected), V),
            native_active_strand_equivalent_density_Ainv3=ratio(
                sum(eq["active"][i] for i in selected), V
            ),
            native_active_fraction_of_strand_equivalents=ratio(
                sum(eq["active"][i] for i in selected),
                sum(eq["all"][i] for i in selected),
            ),
            two_core_candidate_fraction=ratio(
                sum(eq["core"][i] for i in selected),
                sum(eq["all"][i] for i in selected),
            ),
            dangling_fraction=ratio(
                sum(eq["dangling"][i] for i in selected),
                sum(eq["all"][i] for i in selected),
            ),
            self_loop_fraction=ratio(
                sum(eq["self_loop"][i] for i in selected),
                sum(eq["all"][i] for i in selected),
            ),
            isolated_fraction=ratio(
                sum(eq["isolated"][i] for i in selected),
                sum(eq["all"][i] for i in selected),
            ),
            mean_active_junction_degree=ratio(
                sum(degree_sum[i] for i in selected),
                sum(node_count[i] for i in selected),
            ),
        )
        for name, (possible, taken) in reactive.items():
            summary[region][name + "_site_conversion"] = ratio(
                sum(taken[i] for i in selected), sum(possible[i] for i in selected)
            )
        for name, values in component_counts.items():
            summary[region][name + "_fraction_of_matrix_beads"] = ratio(
                sum(values[i] for i in selected), sum(count[i] for i in selected)
            )
    summary["interface_over_reference"] = {
        k: ratio(v, summary["reference"].get(k))
        for k, v in summary["interface"].items()
        if "density" in k or "fraction" in k
    }
    paired = []
    for p in range(len(centers)):
        local = [r for r in per if r["np_molecule"] == p + info["np_molecule_id_start"]]
        values = {}
        for region in ("interface", "reference"):
            chosen = [r for r in local if r["region"] == region]
            values[region] = ratio(
                sum(r["native_active_strand_equivalents"] for r in chosen),
                sum(r["all_strand_equivalents"] for r in chosen),
            )
        if all(v is not None for v in values.values()):
            paired.append(
                dict(
                    np_molecule=p + info["np_molecule_id_start"],
                    **values,
                    difference=values["interface"] - values["reference"],
                )
            )
    write_tsv(out / "paired_particle_active_fractions.tsv", paired)
    if paired:
        differences = np.array([r["difference"] for r in paired])
        summary["particle_spatial_variation"] = dict(
            paired_particles=len(paired),
            median_active_fraction_difference=float(np.median(differences)),
            difference_q25_q75=np.quantile(differences, [0.25, 0.75]).tolist(),
            particles_with_lower_interface_active_fraction=sum(
                r["difference"] < 0 for r in paired
            ),
            interpretation="Spatial variation within one periodic realization; particles are not independent experimental replicas.",
        )
    sensitivity = []
    for lo, hi in c["comparison_windows_A"]:
        for ref in c["comparison_reference_gaps_A"]:
            if ref < hi or any(v not in c["shell_edges_A"] for v in (lo, hi, ref)):
                continue
            near = [
                i
                for i, r in enumerate(rows)
                if r["gap_lower_A"] is not None and lo <= r["gap_lower_A"] < hi
            ]
            far = [
                i
                for i, r in enumerate(rows)
                if r["gap_lower_A"] is not None and r["gap_lower_A"] >= ref
            ]
            a = ratio(
                sum(eq["active"][i] for i in near), sum(eq["all"][i] for i in near)
            )
            b = ratio(sum(eq["active"][i] for i in far), sum(eq["all"][i] for i in far))
            rho_a = ratio(sum(count[i] for i in near), sum(volume[i] for i in near))
            rho_b = ratio(sum(count[i] for i in far), sum(volume[i] for i in far))
            sensitivity.append(
                dict(
                    interface_lower_A=lo,
                    interface_upper_A=hi,
                    reference_minimum_A=ref,
                    near_active_fraction=a,
                    reference_active_fraction=b,
                    active_fraction_ratio=ratio(a, b),
                    matrix_density_ratio=ratio(rho_a, rho_b),
                )
            )
    write_tsv(out / "shell_sensitivity.tsv", sensitivity)
    summary.update(
        matrix_beads_inside_nominal_NP_spheres=int(count[0]),
        proximity_classes=dict(
            collections.Counter(r["proximity_class"] for r in strand_rows)
        ),
        volume_normalization="Disjoint nearest-NP-surface shells sampled uniformly in the periodic current box; shell overlap and other NP interiors counted once.",
        density_definition="Each effective strand contributes 1/N_contour to each contour bead, so integrating its local strand-equivalent density gives one strand. Shared contour endpoints retain native incidence counting.",
        surface_definition="Distance to reconstructed NP center minus metadata radius of the sphere of bead centers; not a calibrated solvent-accessible or silica surface.",
        reference_definition="Operational farther-from-NP comparison region; not proof of unperturbed bulk.",
        uncertainty="Monte Carlo SE describes volume integration only. Per-particle profiles describe spatial heterogeneity within one realization; no independent-replica confidence interval.",
    )
    pairs = []
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            distance = float(np.linalg.norm(mic(centers[j] - centers[i], s["L"])))
            pairs.append(
                dict(
                    np_molecule_i=i + info["np_molecule_id_start"],
                    np_molecule_j=j + info["np_molecule_id_start"],
                    center_distance_A=distance,
                    nominal_surface_gap_A=distance
                    - 2 * info["np_radius_bead_centers_A"],
                )
            )
    write_tsv(out / "NP_pair_gaps.tsv", pairs)
    summary["minimum_NP_center_surface_gap_A"] = min(
        (r["nominal_surface_gap_A"] for r in pairs), default=None
    )
    summary["strand_proximity_by_native_status"] = {
        state: dict(
            collections.Counter(
                r["proximity_class"] for r in strand_rows if r["status"] == state
            )
        )
        for state in sorted(set(r["status"] for r in strand_rows))
    }
    write_json(out / "summary.json", summary)
    write_json(
        out / "np_centers.json",
        dict(centers_A=centers.tolist(), box_lengths_A=s["L"].tolist()),
    )
    return summary, rows


def topology_fingerprint(s):
    import hashlib

    return hashlib.sha256(
        json.dumps(
            sorted((t, min(i, j), max(i, j)) for _, t, i, j in s["bonds"])
        ).encode()
    ).hexdigest()


def dump_frames(path, source, fixed_box=True):
    expected = len(source["pos"])
    expected_ids = np.arange(1, expected + 1)
    with Path(path).open() as f:
        previous = None
        while True:
            first = f.readline()
            if not first:
                break
            if first.strip() != "ITEM: TIMESTEP":
                raise ValueError("Incomplete trajectory header")
            step = int(f.readline())
            if previous is not None and step <= previous:
                raise ValueError("Non-increasing trajectory clock")
            previous = step
            if (
                f.readline().strip() != "ITEM: NUMBER OF ATOMS"
                or int(f.readline()) != expected
            ):
                raise ValueError("Trajectory atom count mismatch")
            if f.readline().strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError(
                    "Rigid-NP interface dynamics requires orthorhombic periodic bulk"
                )
            bounds = np.array(
                [list(map(float, f.readline().split())) for _ in range(3)]
            )
            L = bounds[:, 1] - bounds[:, 0]
            if fixed_box and not np.allclose(L, source["L"], rtol=1e-7, atol=1e-5):
                raise ValueError(
                    "Dynamics trajectory box changed; fixed-volume MSD required"
                )
            columns = f.readline().split()[2:]
            if not set(("id", "mol", "type", "x", "y", "z", "ix", "iy", "iz")).issubset(
                columns
            ):
                raise ValueError("Missing trajectory identity/coordinates/image flags")
            lines = [f.readline() for _ in range(expected)]
            if any(len(line.split()) != len(columns) for line in lines):
                raise ValueError("Truncated or malformed trajectory atom row")
            rows = np.fromstring("".join(lines), sep=" ")
            if rows.size != expected * len(columns):
                raise ValueError("Incomplete trajectory numeric record")
            rows = rows.reshape(expected, len(columns))
            ids = rows[:, columns.index("id")]
            if not np.array_equal(ids, expected_ids):
                raise ValueError("Dump atom IDs must be sorted and complete")
            if not np.array_equal(
                rows[:, columns.index("mol")], source["mol"]
            ) or not np.array_equal(rows[:, columns.index("type")], source["typ"]):
                raise ValueError("Trajectory molecule/type identity changed")
            pos = rows[:, [columns.index(a) for a in "xyz"]]
            images = rows[:, [columns.index("i" + a) for a in "xyz"]]
            if (
                not np.isfinite(pos).all()
                or not np.equal(images, np.rint(images)).all()
            ):
                raise ValueError("Nonfinite coordinate or invalid image flag")
            yield step, pos, images, L, bounds[:, 0]


def trajectory_analysis(tr, s, info, weights, masks, c, out):
    if not Path(tr["path"]).is_file():
        return dict(status="missing_output", path=tr["path"])
    out.mkdir()
    rng = np.random.default_rng(c["seed"])
    matrix_ids = np.flatnonzero(s["typ"] != info["np_atom_type"])
    # A uniform sample is used for shell averages; special midpoint/attachment
    # cohorts are added separately and never mixed into the uniform population.
    uniform = np.sort(
        rng.choice(
            matrix_ids, min(c["max_dynamic_atoms"], len(matrix_ids)), replace=False
        )
    )
    extra = np.flatnonzero(masks["strand_midpoint"] | masks["junction_attachment"])
    selected = np.union1d(uniform, extra)
    uniform_mask = np.isin(selected, uniform)
    group_masks = {"matrix": uniform_mask}
    for name in ("native_active", "two_core_candidate", "defective_strand"):
        group_masks[name] = masks[name][selected] & uniform_mask
    for name in ("strand_midpoint", "junction_attachment"):
        group_masks[name] = masks[name][selected]
    np.save(out / "selected_atom_ids.npy", selected + 1)
    positions = []
    centers_series = []
    com_series = []
    steps = []
    last = None
    last_images = None
    last_centers = None
    unwrapped = None
    unwrapped_centers = None
    com = np.zeros(3)
    flag_mismatches = 0
    matrix_flag_mismatches = 0
    maxjump = 0.0
    rigid_error = 0.0
    reference = None
    mass = s["mass"]
    total_mass = mass.sum()
    L = s["L"]
    for step, pos, images, box, low in dump_frames(tr["path"], s):
        centers, body_reference = body_centers(pos, s["mol"], info, L)
        if reference is None:
            reference = body_reference
        else:
            rigid_error = max(
                rigid_error,
                max(
                    float(np.abs(a - b).max())
                    for a, b in zip(reference, body_reference)
                ),
            )
        if rigid_error > c["rigidity_tolerance_A"]:
            raise ValueError("NP body shape/image reconstruction changed in trajectory")
        if last is None:
            unwrapped = pos[selected].copy()
            unwrapped_centers = centers.copy()
        else:
            delta = mic(pos - last, L)
            maxjump = max(maxjump, float(np.max(np.abs(delta))))
            if maxjump > c["max_frame_displacement_A"]:
                raise ValueError(
                    "Frame interval too coarse for the configured temporal unwrapping bound"
                )
            image_delta = pos - last + (images - last_images) * L
            mismatch = np.any(np.abs(image_delta - delta) > 1e-3, axis=1)
            flag_mismatches += int(mismatch.sum())
            matrix_flag_mismatches += int(mismatch[matrix_ids].sum())
            unwrapped += delta[selected]
            unwrapped_centers += mic(centers - last_centers, L)
            com += (delta * mass[:, None]).sum(axis=0) / total_mass
        positions.append(unwrapped.astype(np.float32).copy())
        centers_series.append(unwrapped_centers.copy())
        com_series.append(com.copy())
        steps.append(step)
        last = pos
        last_images = images
        last_centers = centers
    if len(steps) < 2:
        raise ValueError("At least two complete trajectory frames required")
    stride = np.diff(steps)
    if not (stride == stride[0]).all():
        raise ValueError("Incomplete or irregular trajectory sampling")
    timestep = tr.get("timestep_fs", c["timestep_fs"])
    interval_ps = float(stride[0] * timestep / 1000)
    if timestep <= 0 or not math.isfinite(timestep):
        raise ValueError("Invalid trajectory timestep")
    xyz = np.array(positions, dtype=np.float64)
    centers = np.array(centers_series)
    com = np.array(com_series)
    nframes = len(xyz)
    rows = []
    missing = []
    laglist = c["dw_lags_ps" if tr["kind"] == "dw" else "msd_lags_ps"]
    nb = len(c["shell_edges_A"]) + 1
    for lag in laglist:
        offset = int(round(lag / interval_ps))
        if (
            offset < 1
            or offset >= nframes
            or not math.isclose(offset * interval_ps, lag, rel_tol=0, abs_tol=1e-7)
        ):
            missing.append(lag)
            continue
        possible = np.arange(0, nframes - offset, c["origin_stride"])
        origins = possible[
            np.unique(
                np.linspace(
                    0, len(possible) - 1, min(len(possible), c["max_origins"])
                ).astype(int)
            )
        ]
        records = {name: np.zeros((nb, 9)) for name in group_masks}
        for origin in origins:
            target = origin + offset
            gap, which, normal = nearest(
                xyz[origin], centers[origin], L, info["np_radius_bead_centers_A"]
            )
            bins = shell_index(gap, c)
            dr = xyz[target] - xyz[origin]
            global_dr = dr - (com[target] - com[origin])
            relative = dr - (centers[target] - centers[origin])[which]
            radial = np.einsum("ij,ij->i", relative, normal) ** 2
            r2 = np.einsum("ij,ij->i", relative, relative)
            global_r2 = np.einsum("ij,ij->i", global_dr, global_dr)
            original_np_gap = (
                np.linalg.norm(mic(xyz[target] - centers[target, which], L), axis=1)
                - info["np_radius_bead_centers_A"]
            )
            contact = gap <= c["interface_gap_A"]
            retained = contact & (original_np_gap <= c["interface_gap_A"])
            data = np.array(
                [
                    np.ones(len(selected)),
                    global_r2,
                    r2,
                    radial,
                    np.maximum(0, r2 - radial),
                    r2 * r2,
                    contact,
                    retained,
                    global_r2 * global_r2,
                ]
            ).T
            for name, chosen in group_masks.items():
                for col in range(data.shape[1]):
                    records[name][:, col] += np.bincount(
                        bins[chosen], weights=data[chosen, col], minlength=nb
                    )
        for name, a in records.items():
            for i in range(nb):
                count = a[i, 0]
                rows.append(
                    dict(
                        trajectory_kind=tr["kind"],
                        group=name,
                        lag_ps=lag,
                        **shell_record(i, c),
                        origins=len(origins),
                        bead_origin_pairs=int(count),
                        MSD_global_COM_corrected_A2=ratio(a[i, 1], count),
                        MSD_relative_to_origin_NP_A2=ratio(a[i, 2], count),
                        MSD_radial_relative_A2=ratio(a[i, 3], count),
                        MSD_tangential_relative_A2=ratio(a[i, 4], count),
                        relative_nonGaussian_parameter=(
                            3 * a[i, 5] * count / (5 * a[i, 2] ** 2) - 1
                            if count and a[i, 2] > 0
                            else None
                        ),
                        endpoint_contact_retention=ratio(a[i, 7], a[i, 6]),
                    )
                )
    write_tsv(out / "radial_dynamics.tsv", rows)
    np_msd = []
    for lag in laglist:
        offset = int(round(lag / interval_ps))
        if (
            offset < 1
            or offset >= nframes
            or not math.isclose(offset * interval_ps, lag, abs_tol=1e-7)
        ):
            continue
        origins = np.unique(
            np.linspace(
                0, nframes - offset - 1, min(nframes - offset, c["max_origins"])
            ).astype(int)
        )
        displacement = (
            centers[origins + offset]
            - centers[origins]
            - (com[origins + offset] - com[origins])[:, None, :]
        )
        np_msd.append(
            dict(
                lag_ps=lag,
                origins=len(origins),
                NPs=info["np_count"],
                NP_COM_MSD_A2=float(np.mean(np.sum(displacement**2, axis=2))),
            )
        )
    write_tsv(out / "NP_COM_MSD.tsv", np_msd)
    result = dict(
        status="complete",
        source=tr["path"],
        source_sha256=sha(tr["path"]),
        frames=nframes,
        first_step=steps[0],
        last_step=steps[-1],
        duration_ps=(steps[-1] - steps[0]) * timestep / 1000,
        sample_interval_ps=interval_ps,
        requested_lags_ps=laglist,
        unavailable_lags_ps=missing,
        uniform_sample_atoms=len(uniform),
        special_cohort_atoms=len(selected) - len(uniform),
        maximum_per_frame_component_displacement_A=maxjump,
        image_flag_increment_mismatch_observations=flag_mismatches,
        matrix_image_flag_increment_mismatch_observations=matrix_flag_mismatches,
        maximum_NP_distance_change_A=rigid_error,
        lag_status={
            str(lag): "insufficient_sampling" if lag in missing else "complete"
            for lag in laglist
        },
        coordinate_method="Temporal minimum-image increments of wrapped coordinates; bounded per-frame displacements checked for all atoms. Inherited image flags compared and reported, not used to assemble a network COM.",
        shell_method="Lagrangian cohort: shell and nearest NP assigned at each time origin; relative motion subtracts that same NP COM displacement, without subtracting NP rotation.",
        retention_method="Endpoint proximity to the original NP, conditional on proximity at the origin; not continuous residence or an adsorption lifetime.",
        quality="Network MSD/localization and heterogeneity diagnostics; no asymptotic diffusion or Debye-Waller plateau asserted. Non-Gaussian values can include interfacial anisotropy. Overlapping time origins are correlated.",
    )
    write_json(out / "sampling.json", result)
    return result


def tensile_interface_series(path, axis, source, info, weights, c, out):
    """Eulerian shell occupancy under deformation, with frozen covalent labels."""
    if not Path(path).is_file():
        return dict(status="missing_output", source=str(path))
    out.mkdir()
    matrix = source["typ"] != info["np_atom_type"]
    _, reference = body_centers(source["pos"], source["mol"], info, source["L"])
    rows = []
    sampled = []
    next_strain = 0.0
    max_error = 0.0
    total_frames = 0
    last = None
    k = "xyz".index(axis)

    def analyze(step, pos, L, low):
        current = {**source, "pos": pos, "L": L, "low": low}
        centers, _ = body_centers(pos, source["mol"], info, L)
        gap, _, _ = nearest(pos, centers, L, info["np_radius_bead_centers_A"])
        bins = shell_index(gap, c)
        volume, se, _ = volumes(current, centers, info, c)
        affine = low + (source["pos"] - source["low"]) / source["L"] * L
        residual = mic(pos - affine, L)
        residual -= np.average(residual[matrix], axis=0, weights=source["mass"][matrix])
        residual2 = np.einsum("ij,ij->i", residual, residual)
        hist = lambda a: np.bincount(bins, weights=a, minlength=len(volume))
        count = hist(matrix)
        all_strands = hist(weights["all"])
        active = hist(weights["active"])
        core = hist(weights["core"])
        dangling = hist(weights["dangling"])
        displacement = hist(residual2 * matrix)
        strain = float(L[k] / source["L"][k] - 1)
        sampled.append(dict(step=step, engineering_strain=strain))
        for i, V in enumerate(volume):
            rows.append(
                dict(
                    axis=axis,
                    step=step,
                    engineering_strain=strain,
                    **shell_record(i, c),
                    volume_A3=float(V),
                    volume_MC_SE_A3=float(se[i]) if math.isfinite(se[i]) else None,
                    matrix_beads=int(count[i]),
                    matrix_number_density_Ainv3=ratio(count[i], V),
                    native_active_strand_equivalent_density_Ainv3=ratio(active[i], V),
                    native_active_fraction_of_strand_equivalents=ratio(
                        active[i], all_strands[i]
                    ),
                    two_core_candidate_fraction=ratio(core[i], all_strands[i]),
                    dangling_fraction=ratio(dangling[i], all_strands[i]),
                    minimum_image_affine_residual_MSD_A2=ratio(
                        displacement[i], count[i]
                    ),
                )
            )

    for step, pos, images, L, low in dump_frames(path, source, fixed_box=False):
        total_frames += 1
        _, body = body_centers(pos, source["mol"], info, L)
        max_error = max(
            max_error,
            max(float(np.max(np.abs(a - b))) for a, b in zip(body, reference)),
        )
        if max_error > c["rigidity_tolerance_A"]:
            raise ValueError("NP shape changed during tensile trajectory")
        strain = float(L[k] / source["L"][k] - 1)
        if strain >= next_strain - 1e-7:
            analyze(step, pos, L, low)
            next_strain = (
                math.floor((strain + 1e-7) / c["tensile_profile_interval_strain"]) + 1
            ) * c["tensile_profile_interval_strain"]
        last = (step, pos, L, low)
    if last is None:
        raise ValueError("Empty tensile trajectory")
    if not sampled or sampled[-1]["step"] != last[0]:
        analyze(*last)
    write_tsv(out / "deformation_profiles.tsv", rows)
    result = dict(
        status="complete",
        source=str(path),
        source_sha256=sha(path),
        frames_checked=total_frames,
        sampled_frames=sampled,
        maximum_NP_distance_change_A=max_error,
        membership="Current shell membership (Eulerian); native covalent strand labels frozen from the undeformed snapshot. Redistribution is not creation or destruction of chemical strands.",
        displacement="Minimum-image residual from the homogeneous diagonal cell deformation, with matrix mean translation removed. Includes thermal motion; not Falk-Langer D2min, a local modulus, or a multi-frame diffusion measurement.",
    )
    write_json(out / "sampling.json", result)
    return result


def part2(c):
    out = Path(c["output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    info = json.loads(Path(c["info"]).read_text())
    driver = static_driver()
    normalized = driver.normalize(info)
    if normalized["source_system"] != "composite":
        raise ValueError("NP interface mode requires rigid-composite metadata")
    result = dict(
        status="running",
        source={
            "data": c["data"],
            "data_sha256": sha(c["data"]),
            "info": c["info"],
            "info_sha256": sha(c["info"]),
        },
        config=c,
        analyses=[],
        software={
            "interface_module_sha256": sha(__file__),
            "numpy_version": np.__version__,
            "contour_backend_sha256": sha(
                backend("siliconelab_interface_network", c.get("backend_dir"))
            ),
        },
    )
    write_json(out / "analysis.json", result)
    write_json(out / "resolved_config.json", c)
    try:
        static = out / "static"
        command = [
            sys.executable,
            str(Path(driver.__file__).resolve()),
            c["data"],
            c["info"],
            "--output-dir",
            str(static),
            "--z1-mode",
            "export",
        ]
        if c.get("backend_dir"):
            command += ["--backend-dir", c["backend_dir"]]
        invoke(command, out)
        mapping = out / "network_mapping"
        invoke(
            [
                backend("siliconelab_interface_network", c.get("backend_dir")),
                str(static / "matrix.data"),
                str(static / "normalized_info.json"),
                str(mapping),
            ],
            out,
        )
        s = snapshot(c["data"], info["total_beads"])
        contours, nodes, weights, masks, network = network_map(mapping, s, normalized)
        write_json(out / "network_definitions.json", network)
        result["equilibrium"], _ = profile_snapshot(
            s, info, normalized, contours, nodes, weights, c, out / "equilibrium"
        )
        result["network"] = network
        result["analyses"].append(
            dict(
                analysis="NP-distance network structure, conversion and defects",
                status="complete",
            )
        )
        result["tensile_endpoints"] = {}
        for axis, path in c["tensile_data"].items():
            deformed = snapshot(path, info["total_beads"])
            if (
                topology_fingerprint(s) != topology_fingerprint(deformed)
                or not np.array_equal(s["mol"], deformed["mol"])
                or not np.array_equal(s["typ"], deformed["typ"])
            ):
                raise ValueError(
                    "Covalent topology/identity changed during tensile loading"
                )
            # Reuse native bond-by-bond junction reconstruction in the deformed cell.
            from siliconelab_composite_analysis import project_matrix

            projection = out / ("tensile_" + axis + ".matrix.data")
            project_matrix(path, projection, info, normalized)
            current_mapping = out / ("tensile_" + axis + ".network_mapping")
            invoke(
                [
                    backend("siliconelab_interface_network", c.get("backend_dir")),
                    str(projection),
                    str(static / "normalized_info.json"),
                    str(current_mapping),
                ],
                out,
            )
            current_nodes = read_tsv(current_mapping / "nodes.tsv")
            result["tensile_endpoints"][axis], _ = profile_snapshot(
                deformed,
                info,
                normalized,
                contours,
                current_nodes,
                weights,
                c,
                out / ("tensile_" + axis),
            )
            result["tensile_endpoints"][axis]["source_sha256"] = sha(path)
        result["tensile_evolution"] = {
            axis: tensile_interface_series(
                path, axis, s, info, weights, c, out / ("deformation_" + axis)
            )
            for axis, path in c["tensile_trajectories"].items()
        }
        result["dynamics"] = {}
        for tr in c["trajectories"]:
            label = tr.get("label", tr["kind"])
            if "/" in label or label in result["dynamics"]:
                raise ValueError("Unique simple trajectory labels required")
            result["dynamics"][label] = trajectory_analysis(
                tr, s, info, weights, masks, c, out / ("dynamics_" + label)
            )
        result["analyses"].extend(
            [
                dict(
                    analysis="NP-obstacle primitive paths / entanglements",
                    status="not_implemented",
                    reason="Matrix-only Z1 does not encode rigid particle obstacles.",
                ),
                dict(
                    analysis="Planar film surfaces and surface stress",
                    status="not_applicable",
                    reason="This case is a fully periodic bulk composite.",
                ),
            ]
        )
        if (
            sha(c["data"]) != result["source"]["data_sha256"]
            or sha(c["info"]) != result["source"]["info_sha256"]
        ):
            raise ValueError("Source changed during analysis")
        result["status"] = "complete"
        write_structure_report(result, out)
        if c["plot"]:
            result["plots"] = plot_profiles(out)
    except Exception as e:
        result.update(status="failed", error=str(e))
        raise
    finally:
        write_json(out / "analysis.json", result)
    return result


def write_structure_report(result, out):
    summary = result["equilibrium"]
    lines = [
        "# NP/matrix interface structure",
        "",
        "Source snapshot: `" + result["source"]["data"] + "`",
        "",
        "Distances refer to the metadata sphere of NP bead centers. The reference region is operational, not a separate bulk simulation.",
        "",
        "| Quantity | Interface | Reference |",
        "|---|---:|---:|",
    ]
    for key in (
        "matrix_number_density_Ainv3",
        "native_active_strand_equivalent_density_Ainv3",
        "native_active_fraction_of_strand_equivalents",
        "dangling_fraction",
        "self_loop_fraction",
        "crosslinkers_fraction_of_matrix_beads",
        "strands_site_conversion",
    ):
        values = [summary[region][key] for region in ("interface", "reference")]
        lines.append(
            "| "
            + key
            + " | "
            + " | ".join(f"{v:.6g}" if v is not None else "unavailable" for v in values)
            + " |"
        )
    lines += [
        "",
        "## Graph definitions and quality",
        "",
        "Native strand counts: `"
        + json.dumps(result["network"]["native_status_counts"])
        + "`.",
        "",
        result["network"]["caution"],
        "",
        "Periodicity is checked with corrected cycle winding; the winding rank is "
        + str(result["network"]["corrected_periodic_cycle_winding_rank"])
        + ". Native provisional directional shortest-path tables remain preserved.",
        "",
        summary["density_definition"],
        "",
        summary["uncertainty"],
        "",
        "See `equilibrium/shell_sensitivity.tsv` and `equilibrium/paired_particle_active_fractions.tsv` for boundary and spatial sensitivity.",
        "",
        "## Dynamics",
        "",
    ]
    if not result["dynamics"]:
        lines.append(
            "No trajectory was requested in this static report. This report does not claim completed long-time mobility or tensile analysis."
        )
    for label, dynamics in result["dynamics"].items():
        lines.append(
            "- "
            + label
            + ": "
            + dynamics["status"]
            + (
                "; unavailable requested lags (ps): "
                + str(dynamics["unavailable_lags_ps"])
                if dynamics.get("unavailable_lags_ps")
                else ""
            )
        )
    lines += [
        "",
        "NP-obstacle primitive paths/entanglements remain unimplemented. Noncovalent NP contacts do not create chemical network junctions.",
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))


def tensile_summary(path, axis, c, out):
    raw = np.loadtxt(path, comments="#", ndmin=2)
    if (
        raw.shape[1] != 16
        or len(raw) < 3
        or not np.isfinite(raw).all()
        or (np.diff(raw[:, 0]) <= 0).any()
        or (np.diff(raw[:, 1]) <= 0).any()
    ):
        raise ValueError("Incomplete or nonfinite 16-column bulk tensile table")
    k = "xyz".index(axis)
    strain = raw[:, 4 + k]
    stress = raw[:, 15]
    if (np.diff(strain) <= 0).any():
        raise ValueError("Axial strain must increase")
    if not np.allclose(raw[:, 14], -0.101325 * raw[:, 7 + k], rtol=1e-5, atol=1e-5):
        raise ValueError("Tensile stress sign/unit contract failed")
    bins = np.floor(strain / c["tensile_strain_block"]).astype(int)
    curves = []
    for b in np.unique(bins):
        selected = bins == b
        curves.append(
            dict(
                strain=float(strain[selected].mean()),
                engineering_stress_MPa=float(stress[selected].mean()),
                stress_sample_SD_MPa=(
                    float(stress[selected].std(ddof=1)) if selected.sum() > 1 else None
                ),
                samples=int(selected.sum()),
            )
        )
    write_tsv(out / ("tensile_" + axis + "_curve.tsv"), curves)
    fits = []
    for lo, hi in c["tensile_fit_windows"]:
        selected = [r for r in curves if lo <= r["strain"] <= hi]
        if len(selected) < 3:
            fits.append(dict(window=[lo, hi], status="insufficient_sampling"))
            continue
        x = np.array([r["strain"] for r in selected])
        y = np.array([r["engineering_stress_MPa"] for r in selected])
        slope, intercept = np.polyfit(x, y, 1)
        res = y - (slope * x + intercept)
        sst = float(np.sum((y - y.mean()) ** 2))
        fits.append(
            dict(
                window=[lo, hi],
                status="descriptive_fit",
                Young_modulus_MPa=float(slope),
                intercept_MPa=float(intercept),
                R2=1 - float(res @ res) / sst if sst else None,
                strain_blocks=len(x),
            )
        )
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    result = dict(
        status="complete",
        axis=axis,
        source=str(path),
        source_sha256=sha(path),
        samples=len(raw),
        last_sample_mean_strain=float(strain[-1]),
        last_sample_mean_time_ns=float(raw[-1, 1] / 1e6),
        small_strain_fits=fits,
        tensile_work_density_MJ_m3=float(integrate(stress, strain)),
        stress_definition="Saved bulk engineering stress; tension-positive axial composite virial stress times current transverse area / initial transverse area. No film-thickness correction.",
        quality="Window and block sensitivity are reported; one loading direction in one realization is not an independent replica. OLS slopes do not establish a statistically converged modulus. Final table strain is a block mean.",
    )
    write_json(out / ("tensile_" + axis + "_summary.json"), result)
    return result


def part3(c):
    if not c.get("part2_report"):
        raise ValueError("Composite Part 3 requires a matching Part 2 report")
    p2path = Path(c["part2_report"])
    p2 = json.loads(p2path.read_text())
    if (
        p2.get("status") != "complete"
        or p2["source"]["info_sha256"] != sha(c["info"])
        or p2["source"]["data_sha256"] != sha(c["data"])
    ):
        raise ValueError(
            "Part 2 source metadata/snapshot hash mismatch or incomplete Part 2"
        )
    out = Path(c["output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    result = dict(
        status="running",
        source=p2["source"],
        part2_report={"path": str(p2path), "sha256": sha(p2path)},
        config=c,
        analyses=[],
    )
    write_json(out / "analysis.json", result)
    try:
        result["tensile"] = {
            axis: tensile_summary(path, axis, c, out)
            for axis, path in c["tensile_tables"].items()
        }
        result["interface_structure"] = p2["equilibrium"]
        result["network"] = p2["network"]
        result["interface_tensile_evolution"] = p2.get("tensile_evolution", {})
        coupling = []
        for label, dyn in p2["dynamics"].items():
            if dyn.get("status") != "complete":
                continue
            rows = read_tsv(
                p2path.parent / ("dynamics_" + label) / "radial_dynamics.tsv"
            )
            for lag in sorted(set(float(r["lag_ps"]) for r in rows)):
                r = {}
                for region in ("interface", "reference"):
                    selected = [
                        v
                        for v in rows
                        if v["group"] == "matrix"
                        and v["region"] == region
                        and float(v["lag_ps"]) == lag
                    ]
                    n = sum(int(v["bead_origin_pairs"]) for v in selected)
                    for metric in (
                        "MSD_global_COM_corrected_A2",
                        "MSD_relative_to_origin_NP_A2",
                        "MSD_radial_relative_A2",
                        "MSD_tangential_relative_A2",
                    ):
                        r[region + "_" + metric] = (
                            sum(
                                float(v[metric]) * int(v["bead_origin_pairs"])
                                for v in selected
                                if v[metric]
                            )
                            / n
                            if n
                            else None
                        )
                coupling.append(
                    dict(
                        trajectory=label,
                        lag_ps=lag,
                        **r,
                        relative_MSD_interface_over_reference=(
                            ratio(
                                r["interface_MSD_relative_to_origin_NP_A2"],
                                r["reference_MSD_relative_to_origin_NP_A2"],
                            )
                            if r["interface_MSD_relative_to_origin_NP_A2"] is not None
                            else None
                        ),
                    )
                )
        write_tsv(out / "structure_mobility_coupling.tsv", coupling)
        result["structure_mobility"] = coupling
        result["analyses"] = [
            dict(
                analysis="Bulk independent xyz tensile curves",
                status="complete" if len(result["tensile"]) == 3 else "partial",
            ),
            dict(
                analysis="NP-distance structure/mobility comparison",
                status="complete" if coupling else "missing_output",
            ),
            dict(
                analysis="NP-obstacle primitive paths / entanglements",
                status="not_implemented",
            ),
            dict(
                analysis="Planar surface stress and oil descriptors",
                status="not_applicable",
            ),
        ]
        info = json.loads(Path(c["info"]).read_text())
        interaction = info.get("np_matrix_interaction") or (
            "wca" if "wca" in info.get("np_pair_model", "").lower() else "unspecified"
        )
        result["NP_interactions"] = dict(
            matrix=interaction,
            np_np=info.get("np_np_interaction", interaction),
            lambda_value=info.get("np_lambda"),
        )
        result["interpretation_limits"] = [
            "NP–matrix interaction mode from source metadata: "
            + interaction
            + ". Proximity and displacement observables are not measurements of adhesion free energy.",
            "Contact/proximity to a noncovalent NP does not create a chemical junction or an elastically active strand.",
            "Nearest-surface density, active-strand fraction and mechanical response are distinct observables; this one case cannot establish an interaction-strength trend.",
            "Native-active and giant-component 2-core labels are graph diagnostics. Pendant cycles and NP-induced constraints require additional mechanical or obstacle-aware primitive-path validation.",
            "Reported local density/defects describe this realization; reference shells and correlated time origins are not independent bulk controls or replicas.",
        ]
        result["status"] = "complete"
        if c["plot"]:
            result["plots"] = plot_mechanics(out)
        lines = [
            "# NP/matrix interface analysis",
            "",
            "Source: `" + c["data"] + "`",
            "",
            "| Interface / reference quantity | Ratio |",
            "|---|---:|",
        ]
        for k, v in result["interface_structure"]["interface_over_reference"].items():
            lines.append(
                "| "
                + k
                + " | "
                + (f"{v:.4g}" if v is not None else "unavailable")
                + " |"
            )
        lines += [
            "",
            "## Tensile fits",
            "",
            "| Axis | Fit window | E (MPa) | R² |",
            "|---|---|---:|---:|",
        ]
        for axis, t in result["tensile"].items():
            for fit in t["small_strain_fits"]:
                if "Young_modulus_MPa" in fit:
                    lines.append(
                        f"| {axis} | {fit['window']} | {fit['Young_modulus_MPa']:.4g} | {fit['R2']:.4g} |"
                        if fit["R2"] is not None
                        else f"| {axis} | {fit['window']} | {fit['Young_modulus_MPa']:.4g} | unavailable |"
                    )
        lines += ["", "## Interpretation limits", ""] + [
            "- " + s for s in result["interpretation_limits"]
        ]
        (out / "report.md").write_text("\n".join(lines) + "\n")
    except Exception as e:
        result.update(status="failed", error=str(e))
        raise
    finally:
        write_json(out / "analysis.json", result)
    return result


def plot_profiles(out):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return dict(status="unavailable", reason="matplotlib missing")
    rows = read_tsv(out / "equilibrium/radial_profiles.tsv")
    rows = [
        r
        for r in rows
        if r["gap_lower_A"] and r["gap_upper_A"] and float(r["gap_lower_A"]) >= 0
    ]
    x = np.array(
        [(float(r["gap_lower_A"]) + float(r["gap_upper_A"])) / 2 for r in rows]
    )
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    panels = [
        ("Matrix bead density (Å⁻³)", ["matrix_number_density_Ainv3"]),
        (
            "Strand-equivalent density (Å⁻³)",
            [
                "native_active_strand_equivalent_density_Ainv3",
                "two_core_candidate_strand_equivalent_density_Ainv3",
            ],
        ),
        (
            "Fraction of strand equivalents",
            [
                "native_active_fraction_of_strand_equivalents",
                "two_core_candidate_fraction_of_strand_equivalents",
            ],
        ),
        (
            "Defect fraction of strand equivalents",
            [
                "dangling_fraction_of_strand_equivalents",
                "self_loop_fraction_of_strand_equivalents",
                "isolated_fraction_of_strand_equivalents",
            ],
        ),
    ]
    for ax, (ylabel, keys) in zip(axes.flat, panels):
        for key in keys:
            ax.plot(
                x,
                [float(r[key]) if r[key] else np.nan for r in rows],
                "-o",
                label=key.replace("_fraction_of_strand_equivalents", "")
                .replace("_strand_equivalent_density_Ainv3", "")
                .replace("_", " "),
            )
        ax.set(xlabel="Distance from nominal NP surface (Å)", ylabel=ylabel)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    for ext in ("png", "pdf"):
        fig.savefig(out / ("interface_profiles." + ext), dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    extra = [
        (
            "Crosslinker fraction of matrix beads",
            ["crosslinkers_fraction_of_matrix_beads"],
        ),
        (
            "Reactive-site conversion",
            ["strands_site_conversion", "crosslinkers_site_conversion"],
        ),
        ("Active bond orientation: radial P₂", ["native_active_bond_radial_P2"]),
        ("Mean active junction degree", ["mean_active_junction_degree"]),
    ]
    for ax, (ylabel, keys) in zip(axes.flat, extra):
        for key in keys:
            ax.plot(
                x,
                [float(r[key]) if r[key] else np.nan for r in rows],
                "-o",
                label=key.replace("_", " "),
            )
        ax.set(xlabel="Distance from nominal NP surface (Å)", ylabel=ylabel)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(out / ("interface_chemistry_orientation." + ext), dpi=180)
    plt.close(fig)
    return dict(
        status="complete",
        files=[
            "interface_profiles.png",
            "interface_profiles.pdf",
            "interface_chemistry_orientation.png",
            "interface_chemistry_orientation.pdf",
        ],
    )


def plot_mechanics(out):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return dict(status="unavailable", reason="matplotlib missing")
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for axis in "xyz":
        path = out / ("tensile_" + axis + "_curve.tsv")
        if path.exists():
            rows = read_tsv(path)
            ax.plot(
                [float(r["strain"]) for r in rows],
                [float(r["engineering_stress_MPa"]) for r in rows],
                label=axis,
            )
    ax.set(xlabel="Engineering strain", ylabel="Engineering stress (MPa)")
    ax.legend(title="Independent loading")
    ax.grid(alpha=0.25)
    for ext in ("png", "pdf"):
        fig.savefig(out / ("tensile_xyz." + ext), dpi=180)
    plt.close(fig)
    return dict(status="complete", files=["tensile_xyz.png", "tensile_xyz.pdf"])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--part", choices=["2", "3", "all"], default="all")
    a = p.parse_args(argv)
    c = configure(a.config)
    if a.part == "2":
        part2(c)
    elif a.part == "3":
        part3(c)
    else:
        root = Path(c["output_dir"])
        root.mkdir(parents=True, exist_ok=False)
        p2 = {**c, "output_dir": str(root / "part2")}
        part2(p2)
        part3(
            {
                **c,
                "output_dir": str(root / "part3"),
                "part2_report": str(root / "part2/analysis.json"),
            }
        )


if __name__ == "__main__":
    main()
