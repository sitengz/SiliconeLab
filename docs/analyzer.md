# Generalized structural analyzer

The analyzer accepts an equilibrated LAMMPS data snapshot and its original
generator `.info` file. Geometry and network calculations use C++17; a Python
3.9+ standard-library controller selects the input adapter, runs the C++
programs and external Z1+, validates results, and writes provenance.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
./build/siliconelab_analyzer data.CASE.npt_eq CASE.info
```

The default output is a **new** `analysis_CASE` directory next to the snapshot.
An existing output directory is rejected; use `--output-dir` for another run.
Original snapshots, metadata, and simulation outputs are never modified.

On Nova, run substantial analysis and Z1+ on a compute node. The controller
does not submit jobs. It searches for `Z1+` on PATH, then the explicitly
authorized `/work/wxia/siteng/Z1+` installation. Elsewhere use:

```sh
./build/siliconelab_analyzer data.CASE.npt_eq CASE.info \
    --z1-command /path/to/Z1+ --output-dir results/CASE-analysis
```

Z1+ remains an external dependency and is not redistributed. Its working
directory is `analysis_CASE/z1`, separated from simulation and installation
folders because Z1+ removes/replaces its own result files before a run.

## Inputs and selection

Supported metadata:

Snapshots must use orthorhombic bounds and LAMMPS full-style atoms. Triclinic
boxes are rejected explicitly; tilt factors are never silently ignored.

| Input | Metadata contract | Selection |
| --- | --- | --- |
| Standalone oil | `model = standalone PDMS/PMPS silicone oil` | All oil molecules, including chain-length distributions |
| Elastomer | `pdms-elastomer-model-info`, version 3+ | Linear, ring, star, or grafted network, with optional PDMS filler |
| Coating | `V22-model-info` or `V35-model-info`, version 2 | Bifunctional network, crosslinkers, moderators, PDMS/PMPS/copolymer oil |

The controller saves `source.info` and an explicit `normalized_info.json`;
the original metadata is not relabeled or edited. Coating component order is
network/crosslinker/oil/moderator; elastomer order is
network/crosslinker/moderator/filler. The adapter preserves those original
molecule ranges. Atom type 1 alone cannot distinguish network from oil.

Checks include total and component bead counts, molecule sizes, bead masses,
initial bonds plus created reaction bonds, permitted reaction partners,
functional-site conservation, MPS backbone/pendant counts, and oil backbone
length/composition histograms. Coating balanced MPS allocation can produce
two different bead counts at the same repeat length and is accepted explicitly.
Different seeds with identical composition cannot be distinguished solely by
snapshot counts; the user must supply matching metadata. Both input hashes
are recorded and rechecked at the end.
The manifest also records hashes of the controller, C++ executables and Z1+
launcher. The launcher hash does not identify every file in an external Z1+
installation; preserve that installation separately for reproduction.

## Twelve analyses

| # | Analysis | Systems | Primary output |
| --- | --- | --- | --- |
| 1 | Composition and snapshot validation | All | `snapshot_statistics.tsv`, `analysis.json` |
| 2 | Density and instantaneous kinetic temperature | All | `snapshot_statistics.tsv` |
| 3 | Crosslink conversion | Network systems | `chemical_network_statistics.tsv`, native network statistics |
| 4 | Chemical and reduced-network connectivity | Network systems | Chemical statistics, native topology report |
| 5 | Molecular degrees and reduced-junction properties | Network systems | `reactive_molecule_degrees.tsv`, `junction_properties.CASE.tsv` |
| 6 | Active, dangling, loop and isolated-parent defects | Network systems | `effective_strands.CASE.tsv`, topology/basic reports |
| 7 | Cycle rank, parallel edges, winding and periodic-image shortest paths | Network systems | `directional_self_paths.CASE.tsv`, topology report |
| 8 | Contour length and end-to-end dimensions | All | Native strand tables and `oil_chain_properties.tsv` |
| 9 | Gyration tensor, eigenvalues, Rg, shape, orientation and tortuosity | All | Native strand tables and `oil_chain_properties.tsv` |
| 10 | Affine/phantom structural modulus estimates | Network systems | `network_statistics.CASE.tsv` |
| 11 | Mapped Z1 contour export | All with suitable open contours | `z1/config.Z1`, `z1/chain_mapping.json` |
| 12 | Primitive paths and entanglement observables | Selected Z1 contours | `z1/primitive_paths.tsv`, `z1/primitive_path_summary.json` |

The machine-readable `analysis.json` gives a status for every analysis:
`complete`, `not_applicable`, `not_run`, `blocked`, or `pending`. A fatal
validation or execution error sets the overall status to `failed` and records
the diagnostic. Pure oil reports analyses 3–7 and 10 as `not_applicable`,
rather than inventing a zero-modulus chemical network.

Density uses the complete snapshot box volume. Kinetic temperature requires
complete velocities and removes whole-system center-of-mass motion; absent
velocities produce `null`/`nan`, not a guessed thermostat temperature. Film
box density includes vacuum and should not be labeled material density.

Network conversion separately reports strand, crosslinker, moderator and
limiting-site conversion. Chemical connectivity includes all reactive
molecules and excludes oil. Reduced-network connectivity uses the pinned
elastomer strand/junction reduction; it is a different graph. Full original
network tables are retained for comparisons with the upstream analyzer.

Network conformation uses effective strands. Oil Rg/tensors use all beads,
including MPS pendants, with bead-number weighting. Oil Lc/Ree use its ordered
backbone. Covalent paths are unwrapped bond by bond, so a molecule spanning a
periodic box is not reconstructed by displacement from its first bead. MPS
pendant branches are validated as single attachments to MPS backbone beads.

## Z1 definition and validation

For network-containing mixtures the default selection is `network+oil`:

- Retain the upstream `network` export of active, dangling, dangling-loop and
  self-loop paths whose parent has at least one reacted site.
- Preserve the upstream disjoint contour partition for rings, stars and
  grafted networks, including documented endpoint/center exclusions.
- Append all oil backbones so oil participates in the chosen entanglement
  environment. Pendants, crosslinker geometry and reaction bonds are not
  exported as contour beads.

Use `--z1-selection network` to reproduce the network-only environment.
For standalone oil, all oil backbones are exported regardless of this option.
These selections answer different questions and are recorded with every result.
The export is the disconnected-contour Z1 representation; it is not a
crosslink-preserving primitive-path model of the entire chemical graph.

All contours are in physical angstrom coordinates with the full box preserved;
no hidden scaling is applied. Contours with at most two beads are excluded
from the combined input and recorded because Z1+ defines true chains as having
more than two beads. Duplicate source-bead ownership is rejected.

The importer validates chain count, box lengths, node coordinates and
per-chain N/Lpp/Z files against the mapping and shortest paths. It retains
the native 15-column summary, including tube/Ne estimators, and writes
network/oil group summaries. Undefined native values become JSON `null`.
Kink density counts flag occurrences per exported contour; it does not
deduplicate physical constraints shared between contours.

`strand_properties_with_z1.CASE.tsv` adds mapped Lpp to exported network
strands. Unexported strands remain unavailable. For trimmed ring/star/grafted
paths, Lpp belongs to the export contour rather than the original full
effective strand. Native tables are preserved instead of overwritten.

Native Z1 input does not encode nonperiodic film confinement. A film can
receive static analysis and an export, but automatic Z1 execution is blocked
unless the caller supplies `--film-z1-validated` after validating the installed
package's treatment of that geometry. `--geometry film` explicitly overrides
bulk metadata for a derived oil-film snapshot; other systems use their saved
geometry by default.

To generate inputs without running Z1+:

```sh
./build/siliconelab_analyzer data.CASE.npt_eq CASE.info --z1-mode export
```

Exit codes: `0` for completed requested calculations (including intentional
export-only mode), `2` when Z1 is blocked by a missing executable or unvalidated
film geometry, and `1` for invalid inputs or execution/result-validation failure.
The manifest distinguishes export-only `partial` from all applicable analyses
being complete.

## Limits and verification

These are static snapshot measurements. No trajectory or additional physical
MD time is required for the twelve analyses. Z1+ is an additional computation.
Time-averaged structure, uncertainty, reaction kinetics, equilibration checks,
MSD/diffusion and measured tensile response require other inputs/calculations.
Structural moduli are estimates, not measured mechanical response.
They retain the upstream affine/phantom definitions and use the final target
temperature from the info file; the instantaneous velocity temperature is
reported separately.

Local contract tests cover oil/network applicability, known chain geometry
across a periodic boundary, balanced-MPS mixtures, metadata mismatch, triclinic
rejection, exact native table parity for a reacted network with parallel edges,
output preservation and fail-closed external
Z1 integration. On 2026-10-08, Nova validation completed with exit `0:0` using
the installed Z1+ on all four full-size snapshots. Detailed results and software
hashes are kept in machine-specific local profiles.

| Full-size snapshot | Beads | Applicable analyses passed | Selected Z1 contours |
| --- | ---: | ---: | ---: |
| N16 PDMS oil | 100,000 | 6; six network analyses not applicable | 6,250 oil |
| N40 elastomer | 168,000 | 12 | 2,979 network |
| V22 with N30 PDMS oil | 136,020 | 12 | 888 network + 453 oil |
| V35 with N15 copolymer oil | 136,785 | 12 | 303 network + 638 oil |

```sh
ctest --test-dir build --output-on-failure
```
