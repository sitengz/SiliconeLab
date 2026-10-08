#!/usr/bin/env python3
"""Two-input static analysis controller; C++ does the geometry and graph work."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
from datetime import datetime, timezone

ANALYSES = [
    "composition_validation", "density_temperature", "crosslink_conversion",
    "network_connectivity", "junction_structure", "network_defects",
    "graph_topology", "chain_dimensions", "chain_shape", "structural_modulus",
    "z1_export", "primitive_path_entanglements",
]
NETWORK_ONLY = {3, 4, 5, 6, 7, 10}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def normalize(info, geometry=None):
    """Preserve source semantics while mapping component ranges to the C++ reader."""
    fmt = info.get("format")
    if fmt == "pdms-elastomer-model-info":
        if info.get("format_version", 0) < 3:
            raise ValueError("elastomer metadata version 3 or newer is required")
        system = "elastomer"
        result = json.loads(json.dumps(info))
        comps = result["components"]
        strand = result.get("strand", {
            "topology": "linear", "functionality": 2,
            "reactive_bead_sites": [1, comps["strands"]["N"]],
        })
        result["strand"] = strand
        original_order = ["strands", "crosslinkers", "moderators", "filler"]
        ranges = [(comps[k].get("molecule_id_start"), comps[k].get("molecule_id_end"))
                  for k in original_order if comps[k]["M"]]
        expected = sorted(ranges)
        cursor = 1
        for start, end in expected:
            if start != cursor or end < start:
                raise ValueError("invalid, overlapping or noncontiguous molecule ranges")
            cursor = end + 1
        if cursor - 1 != result["composition"]["total_molecules"]:
            raise ValueError("molecule ranges do not cover metadata total")
        for k in original_order:
            c = comps[k]
            if c["M"] and c["molecule_id_end"] - c["molecule_id_start"] + 1 != c["M"]:
                raise ValueError("component molecule range/count mismatch")
        initial_bonds = info["topology"]["bonds"]
        type2 = info["crosslinker"]["type2_sites_total_including_moderators"]
        type3 = info["crosslinker"]["type3_sites_total"]
        lower = upper = [comps[k]["N"] for k in original_order]
    elif fmt in ("V22-model-info", "V35-model-info"):
        if info.get("format_version") != 2:
            raise ValueError("coating metadata version 2 is required")
        system = "coating"
        originals = info["components"]
        names = {"strands": "network_strands", "crosslinkers": "crosslinkers",
                 "filler": "silicone_oil_filler", "moderators": "star_moderators"}
        comps = {}
        cursor = 1
        for k in ("strands", "crosslinkers", "filler", "moderators"):
            c = dict(originals[names[k]])
            c["molecule_id_start"] = cursor if c["M"] else None
            c["molecule_id_end"] = cursor + c["M"] - 1 if c["M"] else None
            cursor += c["M"]
            comps[k] = c
        oil = info["silicone_oil"]
        if comps["filler"]["M"]:
            n = oil["repeat_units_per_chain"] + oil["mps_repeats_per_chain"]
            extra = int(oil.get("chains_with_extra_mps", 0) > 0)
            comps["filler"]["N"] = n
        else:
            n = extra = 0
        order = ["strands", "crosslinkers", "moderators", "filler"]
        lower = [comps[k]["N"] for k in order]
        upper = list(lower)
        upper[3] = n + extra
        result = {
            "case_name": info["case_name"], "geometry": info["geometry"],
            "components": comps, "composition": info["composition"],
            "strand": {"topology": "linear", "functionality": 2,
                       "reactive_bead_sites": [1, comps["strands"]["N"]]},
            "crosslinker": info["crosslinker"],
            "force_field": {"bond_type_map": {"crosslink": 2}},
            "simulation_template": info["simulation_template"],
        }
        for key in ("film_thickness_angstrom",):
            if key in info:
                result[key] = info[key]
        initial_bonds = info["topology"]["bonds"]
        type2 = info["crosslinker"]["type2_sites_total_including_moderators"]
        type3 = info["crosslinker"]["type3_sites_total"]
        if cursor - 1 != info["composition"]["total_molecules"]:
            raise ValueError("coating component molecule counts do not match total")
    elif info.get("model") == "standalone PDMS/PMPS silicone oil":
        system = "oil"
        oil = info["composition"]
        m = oil["chain_count"]
        beads = info["topology_counts"]["atoms"]
        comps = {k: {"N": 0, "M": 0, "beads": 0,
                     "molecule_id_start": None, "molecule_id_end": None}
                 for k in ("strands", "crosslinkers", "moderators", "filler")}
        comps["filler"] = {"N": 0, "M": m, "beads": beads,
                           "molecule_id_start": 1, "molecule_id_end": m}
        lengths = [c["length"] for c in oil["chain_length_distribution"]]
        if sum(c["chains"] for c in oil["chain_length_distribution"]) != m:
            raise ValueError("oil chain-length histogram/count mismatch")
        lower, upper = [0, 0, 0, min(lengths)], [0, 0, 0, 2 * max(lengths)]
        result = {
            "case_name": info["case_name"], "geometry": info["initial_box"]["geometry"],
            "components": comps, "composition": {"total_beads": beads, "total_molecules": m},
            "strand": {"topology": "linear", "functionality": 0, "reactive_bead_sites": []},
            "crosslinker": {"reactive_bead_sites": []},
            # Oil bond type 2 is covalent MPS backbone, not a reaction bond.
            "force_field": {"bond_type_map": {"crosslink": -1}},
            "simulation_template": info["simulation_template"],
        }
        initial_bonds = info["topology_counts"]["bonds"]["total"]
        type2 = type3 = 0
    else:
        raise ValueError("unsupported info format; expected silicone oil, elastomer v3+, or V22/V35 v2")
    result.update(format="siliconelab-analysis-info", format_version=1,
                  source_system=system, source_format=fmt or info.get("model"),
                  network_enabled=int(system != "oil"), initial_bonds=initial_bonds,
                  initial_type2_sites=type2, initial_type3_sites=type3,
                  minimum_molecule_beads=lower, maximum_molecule_beads=upper)
    order = ["strands", "crosslinkers", "moderators", "filler"]
    result["expected_component_beads"] = [result["components"][k]["beads"] for k in order]
    if sum(result["expected_component_beads"]) != result["composition"]["total_beads"]:
        raise ValueError("component bead counts do not match total")
    if system == "oil":
        types=info["atom_types"]
        result["expected_mps_backbone"]=types["4"]["count"]
        result["expected_mps_pendant"]=types["5"]["count"]
    else:
        types=info["force_field"]["atom_types"]
        mps=info.get("silicone_oil",{}).get("total_mps_repeats",0)
        result["expected_mps_backbone"]=mps
        result["expected_mps_pendant"]=mps
    result["expected_type_masses"]={k:v["mass"] for k,v in types.items()}
    if geometry:
        result["geometry"] = geometry
    if result["geometry"] not in ("bulk", "film"):
        raise ValueError("unknown boundary geometry")
    return result


def table(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def metrics(path):
    result = {}
    for row in table(path):
        key = row.get("quantity") or row.get("property") or row.get("key")
        if not key or key in result:
            raise ValueError(f"missing or duplicate metric name in {path}")
        v = float(row["value"])
        result[key] = v if math.isfinite(v) else None
    return result


def run(command, cwd, log):
    with log.open("a") as stream:
        stream.write("COMMAND " + json.dumps([str(x) for x in command]) + "\n")
        stream.flush()
        p = subprocess.run([str(x) for x in command], cwd=cwd, stdout=stream,
                           stderr=subprocess.STDOUT, check=False)
    if p.returncode:
        raise RuntimeError(f"command exited {p.returncode}; see {log}")


def read_z1(path):
    with Path(path).open() as f:
        count = int(f.readline())
        box = [float(x) for x in f.readline().split()]
        lengths = [int(x) for x in f.readline().split()]
        if len(box) != 3 or len(lengths) != count or any(n < 1 for n in lengths):
            raise ValueError("invalid Z1 input header")
        chains = []
        for n in lengths:
            points = []
            for _ in range(n):
                p = [float(x) for x in f.readline().split()]
                if len(p) != 3 or not all(math.isfinite(x) for x in p):
                    raise ValueError("invalid Z1 coordinates")
                points.append(p)
            chains.append(points)
        if f.read().strip():
            raise ValueError("unexpected trailing Z1 input")
    return box, chains


def export_z1(out, case, network, selection):
    """Keep native strand partitioning, and append oil backbones only when requested."""
    chains, mapping = [], []
    stats = metrics(out / "snapshot_statistics.tsv")
    box = [stats["box_L" + a] for a in "xyz"]
    if network:
        native_box, native = read_z1(out / f"config.{case}.Z1")
        if native_box != box and any(abs(a-b) > 1e-8*max(1, abs(a)) for a,b in zip(native_box,box)):
            raise ValueError("native export box differs from snapshot")
        rows = table(out / f"config.{case}.Z1.map.tsv")
        if len(rows) != len(native):
            raise ValueError("native Z1 chain mapping mismatch")
        for points, row in zip(native, rows):
            chains.append(points)
            mapping.append({"component": "network", "molecule": int(row["parent_molecule"]),
                            "effective_strand_id": int(row["effective_strand_id"]),
                            "status": row["status"], "atom_ids": row["z1_atom_ids"]})
    if not network or selection == "network+oil":
        by_molecule = collections.OrderedDict()
        for row in table(out / "oil_backbones.tsv"):
            mol = int(row["molecule"])
            by_molecule.setdefault(mol, []).append(row)
        for mol, rows in by_molecule.items():
            chains.append([[float(r[a + "_A"]) for a in "xyz"] for r in rows])
            mapping.append({"component": "oil", "molecule": mol, "effective_strand_id": None,
                            "status": "uncrosslinked_chain", "atom_ids": ",".join(r["atom"] for r in rows)})
    excluded = []
    kept, kept_map = [], []
    owners = set()
    for points, row in zip(chains, mapping):
        if len(points) <= 2:
            excluded.append({**row, "reason": "Z1+ true chains require more than two beads"})
            continue
        ids = [int(x) for x in row["atom_ids"].split(",") if x]
        if len(ids) != len(points) or len(set(ids)) != len(ids) or owners.intersection(ids):
            raise ValueError("overlapping or mismatched Z1 contour atom mapping")
        owners.update(ids)
        kept.append(points)
        kept_map.append({**row, "z1_chain_id": len(kept), "beads": len(points)})
    folder = out / "z1"
    folder.mkdir()
    config = folder / "config.Z1"
    with config.open("w") as f:
        f.write(f"{len(kept)}\n" + " ".join(format(x, ".15g") for x in box) + "\n")
        f.write(" ".join(str(len(p)) for p in kept) + "\n")
        for points in kept:
            for p in points:
                f.write(" ".join(format(x, ".15g") for x in p) + "\n")
    save_json(folder / "chain_mapping.json", {"selection": selection if network else "all-oil-backbones",
              "coordinate_scale": 1, "box_A": box, "chains": kept_map, "excluded": excluded,
              "oil_pendants_included": False})
    return folder, kept_map, box


def import_z1(folder, mapping, box):
    """Validate installed single-snapshot format and join results by export order."""
    lines = (folder / "Z1+SP.dat").read_text().splitlines()
    lines = [line.split() for line in lines if line.strip()]
    count = int(lines[0][0])
    actual_box = [float(x) for x in lines[1]]
    if count != len(mapping) or len(lines[0]) != 1 or len(actual_box) != 3:
        raise ValueError("Z1+ shortest-path chain/header mismatch")
    if any(abs(a-b) > 1e-5*max(1, abs(b)) for a,b in zip(actual_box,box)):
        raise ValueError("Z1+ result box differs from physical snapshot box")
    cursor = 2
    results = []
    for chain in mapping:
        n = int(lines[cursor][0]); cursor += 1
        if n < 2:
            raise ValueError("invalid primitive-path node count")
        points = [[float(x) for x in line[:5]] for line in lines[cursor:cursor+n]]
        cursor += n
        if len(points) != n or any(len(p) != 5 or not all(math.isfinite(x) for x in p) for p in points):
            raise ValueError("invalid primitive-path coordinates")
        lc = sum(math.dist(a[:3], b[:3]) for a,b in zip(points,points[1:]))
        z = sum(int(p[4] != 0) for p in points)
        segments = [([b[i]-a[i] for i in range(3)]) for a,b in zip(points,points[1:])]
        orientation = [0.5*(3*s[2]*s[2]/sum(v*v for v in s)-1)
                       for s in segments if sum(v*v for v in s)>0]
        results.append({**chain, "Lpp_A": lc, "kinks": z,
                        "primitive_P2_z": statistics.mean(orientation) if orientation else None})
    if cursor != len(lines):
        raise ValueError("unexpected extra Z1+ snapshots/data")
    for filename, key in (("Lpp_values.dat", "Lpp_A"), ("Z_values.dat", "kinks"),
                          ("N_values.dat", "beads")):
        values = [float(x) for x in (folder / filename).read_text().split()]
        if len(values) != len(results) or not all(math.isfinite(x) for x in values):
            raise ValueError(f"invalid {filename}")
        for row, value in zip(results,values):
            if not math.isclose(row[key], value, rel_tol=2e-4, abs_tol=1e-4):
                raise ValueError(f"{filename} disagrees with mapped shortest path")
            # Retain the native result precision after validating the path geometry.
            row[key] = value
    summary = [float(x) for x in (folder / "Z1+summary.dat").read_text().split()]
    if len(summary) != 15 or int(summary[1]) != len(results):
        raise ValueError("invalid single-snapshot Z1+ summary")
    keys = ["timestep", "true_chains", "mean_original_beads", "mean_Ree2_A2",
            "mean_Lpp_A", "mean_kinks", "coil_tube_diameter_A", "coil_tube_step_A",
            "rms_Lpp_A", "Ne_classical_kink", "Ne_modified_kink", "Ne_classical_coil",
            "Ne_modified_coil", "mean_original_bond_A", "bead_number_density_A-3"]
    native_summary = {k: v if math.isfinite(v) else None for k,v in zip(keys,summary)}
    if results:
        with (folder / "primitive_paths.tsv").open("w") as f:
            fields = [k for k in results[0] if k != "atom_ids"]
            writer = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader(); writer.writerows(results)
    volume = math.prod(box)
    groups = {}
    for label in ("all", "network", "oil"):
        selected = [r for r in results if label == "all" or r["component"] == label]
        if selected:
            groups[label] = {"chains": len(selected),
                "mean_Lpp_A": statistics.mean(r["Lpp_A"] for r in selected),
                "mean_kinks": statistics.mean(r["kinks"] for r in selected),
                "contour_length_density_A-2": sum(r["Lpp_A"] for r in selected)/volume,
                "kink_marker_density_A-3": sum(r["kinks"] for r in selected)/volume}
    save_json(folder / "primitive_path_summary.json", {"native": native_summary, "groups": groups,
              "density_definition": "kink flags counted per contour; not deduplicated physical constraints"})
    return {"native": native_summary, "groups": groups}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("equilibrated_data", type=Path)
    parser.add_argument("info", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--backend-dir", type=Path)
    parser.add_argument("--z1-command", default="Z1+")
    parser.add_argument("--z1-mode", choices=("run", "export"), default="run")
    parser.add_argument("--z1-selection", choices=("network", "network+oil"), default="network+oil")
    parser.add_argument("--geometry", choices=("bulk", "film"))
    parser.add_argument("--film-z1-validated", action="store_true",
                        help="explicitly acknowledge validated Z1+ film boundary handling")
    parser.add_argument("--z1-timeout-seconds", type=int, default=3600)
    args = parser.parse_args(argv)
    data, original = args.equilibrated_data.resolve(), args.info.resolve()
    info = json.loads(original.read_text())
    normalized = normalize(info, args.geometry)
    case = normalized["case_name"]
    if not case or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-+." for c in case):
        raise ValueError("unsafe case name")
    out = (args.output_dir or data.parent / f"analysis_{case}").resolve()
    if out.exists():
        raise ValueError(f"output directory already exists: {out}; choose a new directory")
    if not data.is_file():
        raise ValueError("equilibrated data file is missing")
    root = Path(__file__).resolve().parent
    suffix = ".exe" if os.name == "nt" else ""
    default_backend = root if (root / ("siliconelab_static" + suffix)).exists() else root.parent / "build"
    if (root / "Release" / ("siliconelab_static" + suffix)).exists():
        default_backend = root / "Release"
    backend = (args.backend_dir or default_backend).resolve()
    for binary in ("siliconelab_static", "siliconelab_topology"):
        if not (backend / (binary + suffix)).is_file():
            raise ValueError(f"missing C++ backend {backend / (binary + suffix)}; build with CMake first")
    command = shutil.which(args.z1_command)
    if command is None and Path(args.z1_command).is_file():
        command = str(Path(args.z1_command).resolve())
    # Generic installs use PATH; the user's authorized Nova installation is explicit.
    if command is None and args.z1_command == "Z1+" and Path("/work/wxia/siteng/Z1+").is_file():
        command = "/work/wxia/siteng/Z1+"
    out.mkdir(parents=True)
    status = {"schema_version": 1, "case": case, "system": normalized["source_system"],
              "geometry": normalized["geometry"], "started_at": datetime.now(timezone.utc).isoformat(),
              "status": "running", "inputs": {"data": str(data), "info": str(original),
              "data_sha256": sha256(data), "info_sha256": sha256(original)},
              "software": {"controller_sha256": sha256(Path(__file__)),
                  "backends": {name: {"path": str(backend / (name + suffix)),
                      "sha256": sha256(backend / (name + suffix))}
                      for name in ("siliconelab_static", "siliconelab_topology")}},
              "analyses": [{"id": i, "name": name, "status": "pending"}
                           for i,name in enumerate(ANALYSES, 1)], "warnings": []}
    manifest = out / "analysis.json"
    save_json(manifest,status)
    try:
        norm_path = out / "normalized_info.json"
        save_json(norm_path,normalized)
        shutil.copyfile(original,out / "source.info")
        run([backend / ("siliconelab_static" + suffix),data,norm_path,out], out,out / "analysis.log")
        network = bool(normalized["network_enabled"])
        if network:
            run([backend / ("siliconelab_topology" + suffix),data,norm_path,"--output-dir",out],
                out,out / "analysis.log")
        status["snapshot"] = metrics(out / "snapshot_statistics.tsv")
        oil_rows=table(out / "oil_chain_properties.tsv")
        if status["system"] == "oil":
            histogram=collections.Counter(int(r["backbone_beads"]) for r in oil_rows)
            expected={r["length"]:r["chains"] for r in info["composition"]["chain_length_distribution"]}
            if dict(histogram)!=expected:
                raise ValueError("oil backbone length histogram differs from metadata")
        elif status["system"] == "coating" and oil_rows:
            oil=info["silicone_oil"]
            if any(int(r["backbone_beads"])!=oil["repeat_units_per_chain"] for r in oil_rows):
                raise ValueError("coating oil backbone length differs from metadata")
            histogram=collections.Counter(int(r["beads"]) for r in oil_rows)
            base=oil["repeat_units_per_chain"]+oil["mps_repeats_per_chain"]
            extra=oil.get("chains_with_extra_mps",0)
            expected={base:oil["chain_count"]-extra}
            if extra: expected[base+1]=extra
            expected={k:v for k,v in expected.items() if v}
            if dict(histogram)!=expected:
                raise ValueError("balanced MPS allocation differs from metadata")
        for row in status["analyses"][:10]:
            row["status"] = "complete" if network or row["id"] not in NETWORK_ONLY else "not_applicable"
            if row["status"] == "not_applicable":
                row["reason"] = "uncrosslinked oil has no chemical network"
        if status["snapshot"]["velocity_temperature"] is None:
            status["warnings"].append("Temperature unavailable: snapshot has no complete velocities; no value inferred from thermostat setting.")
        if network:
            # Preserve native reduction table and fill Lpp only from a matched Z1 export.
            status["network"] = metrics(out / f"network_statistics.{case}.tsv")
            status["chemical_network"] = metrics(out / "chemical_network_statistics.tsv")
        folder, mapping, box = export_z1(out,case,network,args.z1_selection)
        status["analyses"][10].update(status="complete", chains=len(mapping),
            selection=args.z1_selection if network else "all-oil-backbones")
        status["warnings"].append("Snapshot statistics do not establish equilibration or uncertainty; structural moduli are estimates.")
        if args.z1_mode == "export":
            status["analyses"][11].update(status="not_run",reason="export-only mode requested")
        elif not mapping:
            status["analyses"][11].update(status="not_applicable",reason="no selected contours with more than two beads")
        elif normalized["geometry"] == "film" and not args.film_z1_validated:
            status["analyses"][11].update(status="blocked",reason="native Z1 format does not encode film confinement; validate boundary handling first")
        elif command is None:
            status["analyses"][11].update(status="blocked",reason="Z1+ not found; supply --z1-command")
        else:
            status["z1_command"] = [command,"config.Z1"]
            status["software"]["z1_launcher"] = {"path": command, "sha256": sha256(command)}
            with (folder / "Z1.stdout.log").open("w") as f:
                p = subprocess.run(status["z1_command"],cwd=folder,stdout=f,stderr=subprocess.STDOUT,
                                   timeout=args.z1_timeout_seconds,check=False)
            if p.returncode:
                raise RuntimeError(f"Z1+ exited {p.returncode}; see {folder / 'Z1.stdout.log'}")
            status["primitive_paths"] = import_z1(folder,mapping,box)
            status["analyses"][11]["status"] = "complete"
            if network:
                native = out / f"strand_properties.{case}.tsv"
                rows = table(native)
                values = table(folder / "primitive_paths.tsv")
                by_strand = {r["effective_strand_id"]:r for r in values if r["component"]=="network" and r["effective_strand_id"]}
                for row in rows:
                    matched = by_strand.get(row["strand_id"])
                    if matched:
                        row["Lpp_A"] = matched["Lpp_A"]
                        row["Lpp_source"] = "Z1+ mapped exported contour"
                with (out / f"strand_properties_with_z1.{case}.tsv").open("w") as f:
                    writer=csv.DictWriter(f,list(rows[0]) if rows else [],delimiter="\t")
                    if rows: writer.writeheader();writer.writerows(rows)
                status["warnings"].append("Ring/star/grafted Lpp belongs to the documented trimmed export contour; it is not the original full strand contour.")
        if sha256(data)!=status["inputs"]["data_sha256"] or sha256(original)!=status["inputs"]["info_sha256"]:
            raise RuntimeError("input files changed during analysis")
        blocked = any(r["status"] == "blocked" for r in status["analyses"])
        status["status"] = "partial" if blocked or args.z1_mode=="export" else "complete"
        status["finished_at"] = datetime.now(timezone.utc).isoformat()
        save_json(manifest,status)
        report = [f"SiliconeLab structural analysis: {case}",f"System: {status['system']}; geometry: {status['geometry']}"]
        report += [f"{r['id']:2}. {r['name']}: {r['status']}" + (f" — {r['reason']}" if 'reason' in r else '') for r in status['analyses']]
        report += ["", *status["warnings"]]
        (out / "analysis_report.txt").write_text("\n".join(report)+"\n")
        print(f"Analysis {status['status']}: {manifest}")
        return 2 if blocked else 0
    except Exception as e:
        status.update(status="failed",error=str(e),finished_at=datetime.now(timezone.utc).isoformat())
        save_json(manifest,status)
        raise


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError, KeyError, subprocess.TimeoutExpired) as e:
        print(f"siliconelab_analyzer: {e}",file=sys.stderr)
        sys.exit(1)
