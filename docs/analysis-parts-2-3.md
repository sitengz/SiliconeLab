# Spatial, dynamic, surface and mechanical analysis

Every result preserves the source snapshot, metadata, geometry, selected component,
units, timing, algorithm and sampling limitations. A missing file is distinct from
an inapplicable observable and from insufficient sampling. A successful job proves
execution, not equilibration or scientific convergence.

## Part 2: twelve analyses and their inputs

| Analysis | Applicable populations | Required files / additional sampling |
|---|---|---|
| Component density | All; network/crosslinker/moderator/oil separately | Final data + info; trajectory for averages |
| Surface enrichment | Oil/network mixtures in films | Film data + info, actual material surfaces; trajectory for averages |
| Network profiles | Elastomer/coating | Bonded snapshot + info; local conversion, defects, junctions and strand geometry |
| Surface/core comparisons | Films | Data/trajectory; surface definition, symmetric folding and component selections |
| Phase separation | DMS/MPS mixtures; component mixing separately for PDMS oil | Snapshot/trajectory + repeat/component identities; exclude MPS pendants from repeat counts |
| Concentration structure factor | Distinguishable chemical/component populations | Coordinates + box + labels; report finite-box resolution and shell populations |
| Oil aggregation | Oil populations; distinguish MPS-specific contacts from generic oil contacts | Coordinates + molecule IDs + explicit contact cutoff |
| MSD | Beads by component, mass-weighted oil molecular COM | id/mol/type, wrapped positions + images or unwrapped coordinates, boxes, actual timestep/time mapping |
| Diffusion | Mobile oil | Long trajectory; inspect slope/log-slope, fit range, direction and uncertainty; confined z MSD is not asymptotic diffusion |
| Layer dynamics | All, separately by component | Same trajectory; origin-based layer assignment, film recentering, whole-system COM drift correction |
| Debye–Waller displacement | All, separately by component | Independent 100 ps trajectory every 0.1 ps; plateau selection or explicit lag; inverse displacement is a stiffness proxy |
| Fixed-lag dynamics | All | Exact partner frames; default 10 ns MSD and 10 ps displacement; bead-origin and origin-weighted means are distinct |

A bulk z profile is a uniformity check, not a surface measurement. A pure system has
no composition enrichment. Component identities come from molecule ranges, not just
bead types: PDMS oil and PDMS network may share a bead type. Copolymer oil COM includes
pendant mass; chemical composition counts one marker per backbone repeat.

The completed elastomer/V22/V35 bulk trajectories span 5 ns at 5 ps intervals;
standalone oil has no bulk production trajectory. A 10 ps fixed displacement is
possible from these dumps, while a 10 ns lag and resolved sub-ps Debye–Waller curve
require extra sampling. These original outputs remain preserved.

## Part 3: five analyses and independent property jobs

| Analysis | Applicable system | Required files |
|---|---|---|
| Oil surface tension | Standalone oil film | Equilibration/production pressure tables with time, T, PE, Pxx/Pyy/Pzz, Lx/Ly/Lz and both guard-wall forces; geometry dumps |
| Surface statistics/plots | Same film pressure dataset | Production blocks, drift and wall-contact diagnostics; phase clocks remain separate |
| Mechanical surface stress | Elastomer/coating film | Same pressure/geometry outputs; residual bulk/network stress comparison; not necessarily solid surface free energy |
| Tensile response | Elastomer/coating | Dedicated x/y deformation tables: time, T, engineering strain, Pxx/Pyy/Pzz, cell dimensions and volume; reference box, rate and film material thickness |
| ATSC4i | Oil molecules, standalone or within mixtures | Original generated data/backbone sequence; agreed atomistic end groups; RDKit/Mordred versions; no additional MD |

Surface mechanical estimates use full cell height including vacuum and assume two
surfaces. Nonzero guard-wall contact invalidates an unconfined free-surface result.
Tensile pressure is converted to material stress using a stated material-volume
convention; film vacuum must not dilute the reported stress. The fixed nominal
thickness approximation must be labelled, rather than presented as measured thickness.

## Four matched films and dependency contract

See `workflows/films.json` and `examples/films/*.conf`. The oil film reads a copy of
its completed bulk data and follows the native molecule-preserving conversion.
The three network films are independently generated and cured with matching bulk
compositions/seeds. Cutting a periodic crosslinked network is not used.

The elastomer nominal wall-free thickness is 4 times the measured mean active bulk
strand end-to-end distance (182.84904 Å). Coating thicknesses follow the original
matched-film convention: measured bulk Lz (250.508775126585 and 250.501420989649 Å).
The native coating value is initial full box height, whereas the elastomer value
is nominal material thickness; these conventions remain distinct.

```
verified completed bulk + preparation/preflight
  -> oil film: native 0.5 ns initialization + 50 ns relaxation + 50 ns production
       -> dynamics
  -> each independently generated network film: complete native preparation/cure/cool/MSD
       -> surface: 50 ns relaxation + 50 ns production
            -> dynamics: 5 ns relaxation + 100 ps DW + 25 ns long trajectory
            -> independent tensile-x and tensile-y: 5 ns relaxation + 100 ns loading
```

Surface and tensile jobs start from the same verified parent. Tensile loading uses
engineering strain rate 2e-9 fs^-1 to 20% strain; film normal dimension remains fixed,
and only the transverse lateral direction is pressure controlled. No oil tensile
job is scheduled. Each job starts in its own folder and keeps stdout/stderr there.

Slurm `afterok` and `--kill-on-invalid-dep=yes` enforce scheduler dependencies. A
second gate checks parent completion JSON, final atom counts and required outputs;
changed inputs, missing outputs or LAMMPS errors stop the child. Accepted job IDs
are saved immediately. The submitter checks expanded queues and matching historical
job names before submission; partial graphs resume from the saved state. No automatic
MD retries or shortened production inputs are used. Compute-node preparation performs
force-field preflight; generation, compilation and substantial analysis stay off the
login node. Independent queued/running jobs are limited to eight; the user's explicit
exception allows dependency-linked jobs beyond that count.

Monitoring and analysis distinguish completion from scientific quality. Full stages
must contain their completion marker and required output contracts. The hourly monitor
reports failures, output verification issues and completions, and remains quiet when
there is no meaningful change.

## Two generalized analyzers

Build/install with CMake. `siliconelab_part2` uses the same normalized metadata and
C++ snapshot validation as part 1, a streaming C++ trajectory backend, the pinned
network profile backend, and the pinned coating chemical-phase analyzer when
applicable. NumPy computes operational component composition fields and grid FFTs.
`siliconelab_part3` reads the independent property tables and optional original
chain sequences. Both emit `analysis.json`, detailed tables and a readable `report.md`.

```bash
siliconelab_part2 data.dynamics_final CASE.info --geometry film \
  --trajectory dump.layer_dynamics.lammpstrj \
  --dw-trajectory dump.debye_waller.lammpstrj --output-dir analysis-part2

siliconelab_part3 CASE.info --geometry film \
  --pressure-production pressure.production.dat --pressure-equil pressure.equil.dat \
  --tensile-x ../tensile-x/stress_strain.dat \
  --tensile-y ../tensile-y/stress_strain.dat \
  --material-thickness-A NOMINAL_THICKNESS \
  --initial-data ../data.CASE --geometry-data data.surface_production \
  --part2-report ../analysis-part2/analysis.json --output-dir analysis-part3
```

The oil's native film pressure filenames differ and are mapped by the analysis
runner. Oil has no tensile inputs. The original generated data is retained for
sequence extraction; final reacted atom types must not be used to infer precursor
functional-group sequence.

Part 2 defaults to up to 16 sampled origins separated by 100 trajectory frames.
It pools only exact recorded timestep differences, accepts nonuniform sampling,
uses mass-weighted whole-system COM drift removal, and reports beads separately
from mass-weighted molecular COM. Layer assignment is relative to material COM at
each origin. Diffusion slope fits use 50–90% of the observed duration and remain
preliminary. Global and layer fixed-lag tables use 10 ns and 10 ps; inverse-u2 and
layer/global displacement ratios are included. These origins are correlated, so
an independent-replica uncertainty is not inferred.

Generic component segregation uses mass fraction; DMS/MPS chemistry uses one
backbone marker per repeat. Coarse-grid FFT intensity and the independent-label
noise estimate are explicitly labelled; they are not substituted for the native
chemical concentration-spectrum definition. Both 3D and XY tables are saved.
The lowest-q peak is flagged as box limited. Contact clusters are operational,
with an explicit cutoff, and do not by themselves prove phase separation.
Static fields are single-snapshot measurements; time-averaged structure/coarsening
requires analyzing multiple snapshots and reporting their time windows.

Surface/core regions are obtained from a smoothed total density half-maximum
crossing rather than simulation-cell edges. An unresolved vacuum interval or
multiple disconnected material intervals is reported as insufficient geometry.
Native wall-coordinate network summary tables remain auxiliary after wall release.
Native network conversion, defect and chain geometry profiles retain their own
assignment definitions. Z1 spatial import needs separately validated contour/box
mapping; mixed network/oil primitive-path output is never blindly interpreted as
network-only output.

Part 3 separates phase clocks, checks monotonic times and zero guard-wall forces,
uses nonoverlapping completed blocks, and reports block SEM, an approximate normal
95% interval, and drift. Surface plots require Matplotlib. Tensile reporting fits
stress difference over 0–5% engineering strain and separately reports a
model-dependent reduced-strain fit. Neither is automatically labelled Young's
modulus. Film stress uses a stated nominal material-volume correction.

ATSC4i uses the unchanged native reconstruction with RDKit/Mordred, caches repeated
sequences, and reports each sequence's chain multiplicity. Its provisional end-group
assumption is explicit. Install NumPy/Matplotlib from `requirements-analysis.txt`;
optional descriptor dependencies are RDKit and `mordredcommunity`. Missing packages
are reported, never replaced by invented descriptor values.

Part 3 cross-checks a correlated part-2 report against the exact source metadata
hash. Its surface status also requires resolved vacuum/material interfaces in the
surface final snapshot. A final snapshot check does not establish film integrity
throughout production; the saved geometry trajectory remains required evidence.
