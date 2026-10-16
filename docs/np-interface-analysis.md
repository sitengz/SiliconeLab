# Bulk NP/matrix interface analysis

The static analyzer validates a rigid noncovalent composite and reduces its matrix
covalent graph. The interface extension reuses that exact native reduction, then
maps its contours, junctions and reactive sites to distance from current NP
surfaces. It supports the version-1 generalized and legacy NP metadata formats.
NPs remain excluded from the chemical graph.

Build the C++ targets and install `requirements-analysis.txt`. Use the JSON
configuration [example](../examples/composites/interface-analysis.conf), after
copying it beside the case and setting its paths. Relative paths resolve beside
the configuration file. The output directory must be new.

```sh
python3 src/siliconelab_interface.py --config analysis.conf --part all
```

The combined command creates `OUTPUT/part2` and `OUTPUT/part3`. Existing Part 2
and Part 3 entry points also accept `--config FILE`: each writes directly into
the configured output directory. For separate calls, set a new output directory
for Part 3 and set `part2_report` to Part 2's `analysis.json`. Part 3 requires
matching snapshot and metadata hashes. Legacy positional calls retain their
existing oil/network/film behavior and do not silently run film algorithms on NPs.

## Configurable settings

| Key | Default / units | Meaning |
|---|---|---|
| `data`, `info`, `output_dir` | Required paths | Original undeformed equilibrated composite snapshot, matching full NP metadata, and new output. |
| `backend_dir` | Installed binaries or project build | C++ static, topology and interface-contour backends. |
| `part2_report` | Omitted path | Matching complete Part 2 manifest for a separate Part 3 call. |
| `shell_edges_A` | `[0,5,10,15,20,30,40,60]`, Å | Lower boundaries; also creates an interior bin below zero and an overflow bin above the last edge. |
| `interface_gap_A` | `10`, Å | Main near-interface window, from zero to this boundary. |
| `reference_gap_A` | `30`, Å | Main farther-from-NP comparison region. Both boundaries must be shell edges. |
| `comparison_windows_A` | `[[0,5],[5,10],[0,10],[5,15],[10,20]]`, Å | Sensitivity windows; windows not resolved by the chosen shell edges are omitted. |
| `comparison_reference_gaps_A` | `[20,30,40]`, Å | Reference-distance sensitivity; comparison and reference regions must not overlap. |
| `volume_points`, `volume_replicates` | `250000`, `4` | Uniform random integration points per replicate and number of replicates. |
| `seed` | `20261016` | Reproducible volume integration and uniform atom sampling seed. |
| `trajectories` | Empty list | Entries with `path`, `kind` (`dw`/`msd`), optional unique `label`, and optional `timestep_fs`. |
| `timestep_fs` | `5`, fs | Fallback timestep for trajectory clocks. |
| `dw_lags_ps` | `[0.1,0.5,1,4,10,25,50]`, ps | Explicit short-time displacement lags; no automatic plateau claim. |
| `msd_lags_ps` | `[25,100,500,1000,2500,5000,10000]`, ps | Explicit longer-time displacement lags. Unavailable exact partner frames are reported. |
| `max_dynamic_atoms` | `20000` | Maximum uniform matrix-bead sample. Midpoint and junction-attachment cohorts are added separately. |
| `max_origins`, `origin_stride` | `32`, `1` frames | Cap and eligible spacing of time origins. Eligible origins are sampled across the complete available window. |
| `max_frame_displacement_A` | `20`, Å | Maximum accepted temporal minimum-image displacement component between saved frames. |
| `rigidity_tolerance_A` | `0.02`, Å | Maximum trajectory change of anchor-to-body bead distances, allowing text precision. |
| `tensile_tables`, `tensile_data` | Empty xyz path maps | Bulk 16-column stress tables and final snapshots from independent x/y/z loading. |
| `tensile_trajectories` | Empty xyz path map | Optional complete tensile dumps for interface redistribution and rigidity checks during loading. |
| `tensile_profile_interval_strain` | `0.1` | Engineering-strain interval for saved interface profiles; initial and final frames are included. |
| `tensile_fit_windows` | `[[0.005,0.02],[0.01,0.03],[0.02,0.05]]` | Dimensionless engineering-strain windows for descriptive modulus fits. |
| `tensile_strain_block` | `0.002` | Strain-bin width before fitting/plotting. |
| `plot` | `true` | Write PNG/PDF plots when Matplotlib is available; missing plotting dependency is recorded. |

Resolved parameters, software/input hashes, per-observable applicability and
sampling limitations are saved in `analysis.json`. NumPy is required; plotting
can be disabled. Complete execution does not establish equilibration or convergence.

## Structure, density and graph definitions

NP centers are reconstructed from the current rigid constituent coordinates
using periodic minimum images. Stored initial NP centers are not reused for
later snapshots. Surface gap is center distance minus the metadata radius of
the sphere of bead centers. This is a nominal geometric surface, not a calibrated
solvent-accessible surface. The radius must be small enough for unambiguous local
body reconstruction (diameter below half the shortest cell length).

Every point is assigned to its nearest NP surface. Uniform sampling of the
current periodic box estimates disjoint shell volumes, including overlap with
neighboring shells and other NPs. Volumes sum to the current cell volume and
per-particle volumes sum to pooled volumes. Matrix beads inside nominal NP
spheres are reported separately. Monte Carlo standard errors describe numerical
volume integration, not uncertainty from independent simulation replicas.

`radial_profiles.tsv` reports matrix/component density, junction density and
degree, strand/crosslinker site conversion, native active/dangling/self-loop/
isolated fractions, and native-active bond orientation relative to NP normals.
Physical density and conditional active fraction are separate observables:
steric depletion close to the nominal surface can reduce active-strand density
without an equally large change in network quality.

Each effective strand contributes `1/N_contour` to each of its contour beads.
Its contribution integrates to one strand. Local strand-equivalent densities
therefore recover the native global count when integrated over all shells.
Shared contour endpoints retain native incidence counting. These densities are
not counts of strand midpoints or full strands assigned to their nearest NP.

**Native active** means an effective strand joining two distinct covalent
junctions. The additional **2-core candidate backbone** retains native-active
edges in the largest component after iteratively pruning degree-below-two nodes.
This stricter diagnostic can still retain pendant cycles and does not prove
mechanical activity. Dangling/self-loop labels and actual reaction-site counts
retain their native definitions. Proximity to an NP does not create a covalent
junction. Multiple-NP proximity is an operational geometric classification,
not proof of adsorption, force transmission or an entanglement bridge.

Corrected cycle winding includes the raw node-separation boundary shift and
bond-by-bond junction-to-contour attachments. It is independently tested on a
periodic triangle. It does not retrofit or relabel the native provisional
directional shortest-path tables, which remain preserved in the static output.
No obstacle-aware Z1/primitive-path model has been implemented.

`per_particle_profiles.tsv`, paired active fractions and `shell_sensitivity.tsv`
describe spatial variability and boundary sensitivity. Farther-from-NP shells
are operational reference regions, not guaranteed unaffected bulk controls.

## Dynamics and tensile interpretation

Trajectories must have complete sorted IDs, molecule/type identity, uniform
sampling and a fixed periodic box. Temporal minimum-image coordinate increments
are checked against the configured displacement bound; inherited image-flag
increments are compared and mismatches reported. Sparse frames cannot rule out
multiple unobserved periodic crossings; shorten the sampling interval when this
assumption is questionable. All-atom mass-weighted COM increments remove global
translation. NP centers are reconstructed every frame and NP rigidity is checked.
See the [LAMMPS MSD documentation](https://docs.lammps.org/compute_msd.html) for
the role of COM correction and image flags.

Shell membership and nearest NP are assigned at each time origin. Reported
NP-relative displacement subtracts that same NP's COM motion; it does not remove
NP rotation. Radial/tangential MSD, descriptive non-Gaussian parameters, NP COM
MSD and endpoint proximity retention are reported separately. Endpoint retention
is not continuous residence or an adsorption lifetime. Matrix shell averages
use only the uniform sample; extra strand-midpoint/junction-attachment cohorts
are never mixed into that sample. Origin samples overlap and are correlated.

Network MSD is a mobility/localization measure rather than automatically a
self-diffusion coefficient. Connecting localization and elasticity requires
model-specific assumptions; see [Sorichetti et al., 2023](https://doi.org/10.1063/5.0134271).
The short-time lags are explicit displacement observables; a Debye-Waller plateau
must be assessed rather than assumed from the label `dw`.

Part 3 reads the bulk 16-column engineering/true stress format, including the
composite virial stress and separate matrix/composite temperatures. It applies
no film-thickness correction. Independent xyz curves, strain-block means,
fit-window sensitivity and tensile work are saved. OLS fit quality does not
establish a converged modulus or supply an independent-replica confidence interval.
Deformed endpoint profiles require unchanged covalent topology and identities.
Optional tensile trajectories check complete frame identities and rigid NP shape
throughout loading. Sampled profiles use current shell membership with frozen
undeformed covalent labels; changes describe spatial redistribution, not new
chemical reactions. The minimum-image residual from homogeneous diagonal cell
deformation removes matrix translation and includes thermal motion. It is not
Falk–Langer D2min, a local modulus, or a diffusion measurement.

## One additional dynamics allocation

`scripts/composite_dynamics.py prepare --config dynamics.json` prepares a new
folder; submission is a separate authorized action. Its configuration requires
`source_data`, `source_info`, `forcefield_input` and `output_dir`. The force-field
header is copied from an established cold rigid-composite input through its
first timestep command; no interactions, topology, saved velocities or COM
momentum are changed. Relative paths resolve beside the configuration.

Defaults are 300 K, 5 fs timestep, 50 fs thermostat damping, 100 ps DW sampling
every 0.1 ps and 25 ns MSD sampling every 25 ps. Keys are `temperature_K`,
`timestep_fs`, `thermostat_damping_fs`, `dw_duration_ps`, `dw_sample_ps`,
`msd_duration_ns` and `msd_sample_ps`. The fixed box is inherited from the source;
there is no deformation, barostat, reaction or velocity reset. Integrators stay
active across the two phases. Checkpoints include a named DW boundary restart.

`checkpoint_steps=100000` and `thermo_steps=1000` control output cadence.
Resource defaults (`tasks=96`, `memory=200G`, `time_limit=48:00:00`, `account=wxia`,
`partition=nova`) and the module/logging-helper paths describe the current Nova
profile and can be configured. `independent_job_limit=8`, `lammps=lmp` and
`final_rigidity_tolerance_A=1e-6` are explicit settings. The saved launcher runs a
zero-step compute-node preflight before MD, then verifies complete frame counts,
steps, finite output, source hashes, completion and final NP shape. Existing
evidence is rejected; no automatic MD resubmission is performed.

A completed analysis manifest records execution completion; inspect each observable's
status and `quality_alerts` before using its results. Tensile trajectories are
checked against `rigidity_tolerance_A` before calculating deformation profiles.
If any frame fails, that trajectory has `invalid_rigid_geometry` status, with
per-frame distances in `rigidity_checks.tsv`, offending steps and source hash in
`sampling.json`, and no deformation profiles. Independently validated equilibrium,
final-snapshot, dynamics and bulk stress-table analyses continue. The analyzer
does not repair coordinates or relax the configured tolerance automatically.
