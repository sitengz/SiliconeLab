# Unified bulk composite generation

All four material families use the same executable and config format:

```sh
cmake -S . -B build
cmake --build build --parallel 4
./build/siliconelab_generator --capabilities
./build/siliconelab_generator --help composite
./build/siliconelab_generator --config examples/composites/v22-wca.conf
```

`system = oil`, `elastomer`, and `coating` retain their native settings and
workflows, including films. `system = composite` adds a shared nanoparticle
layer to a selected **periodic bulk** matrix. Composite films are deferred as
requested. The validated NP extension uses the shared generator and static
analyzer. Generation does not execute MD, submit jobs or publish anything.

## Shared file, folder and logging conventions

Composite cases follow the same conventions as oil, elastomer and coating
cases. Use the [existing batch layout](../README.md#shared-software-and-future-runs)
for new work: `runs/<batch>/<case>/configs/model.conf`, `generated/`, applicable
stage folders such as `tensile-x/`, `tensile-y/`, `tensile-z/`, and
`analysis/part1/` for supported static analysis. Batch plans and execution state
remain `workflow.json` and `workflow-state.json`. This is the existing directory
convention, not an additional NP orchestration API.

The generated package keeps `data.<case>`, `in.<case>`, `submit.<case>.sh`,
`<case>.info`, `request.conf`, `generator.conf` and `siliconelab.json`. The NP
metadata/template and preserved `precursor/` files are additions to that common
package. Native simulation output names and checkpoint/stage suffixes are
retained. The static analyzer keeps `analysis.json`, `analysis.log`,
`analysis_report.txt`, `source.info`, `normalized_info.json` and the existing
scientific table names. Select the batch's `analysis/part1/` using `--output-dir`;
without it, the shared CLI default remains `analysis_<case>`.

Use the existing stage launcher and cluster profile for execution and logging.
Generated `submit.<case>.sh` files remain native submission templates, as for
the other systems; generation does not install a Nova-specific logger into them.
Workflow-managed launchers retain their established names (`submit.sbatch` for
simulation stages and `analyze.sbatch` for analysis jobs). Submit from the job's
own folder and keep full scheduler stdout/stderr, LAMMPS output and analyzer
logs there, using the same job-ID naming as the selected launcher. Workflow
stage completion/provenance remains `stage-result.json`.

On Nova, source the existing `nova-logging.sh` helper and use `nova_run` with a
readable purpose and concise summary. Central dated activity logs remain under
`/work/wxia/siteng/log/`; full output stays in the corresponding job folder.
Use saved scripts for routine checks, preserve existing log text, and keep
credentials, encoded transfer payloads and long inline source out of the central
activity log. NP cases do not introduce a second logging scheme.

Historical prototype/test directories are preserved with their original paths,
outputs and job records. Their ad hoc layout is not the convention for future
production batches. Existing oil/elastomer/coating cases are not renamed, and
unsupported composite stages are not enabled by adopting the common layout.

## Choose the matrix

| Composite selection | Matrix settings forwarded unchanged |
|---|---|
| `base = v22` | V22 coating settings: `n1`, `m1`, `n2`, `functionality`, moderators, optional oil recipe, density and seeds. |
| `base = v35` | V35 coating settings, including its folded strand geometry and optional oil recipe. |
| `base = elastomer` | Generic PDMS settings: linear, ring, star or grafted strands; stoichiometry; moderators; neutral filler; conversion-controlled curing. |

Use `base` instead of `model` for composites. The full matrix key tables remain
in [the generator reference](generator.md). These are the currently imported
bases; unsupported names fail explicitly. An additional scientific backend
would require registering its generator and verifying its data/input contract.

Existing `system = coating`, `model = v22|v35` or `system = elastomer` configs
can also opt into the same composite engine by adding NP keys. The original
32-NP configuration still works with a **new** output directory. The previously
inspected data/input under `runs/np-repulsive-v22/v22-np32` remain untouched.

For example:

```ini
system = composite
base = v22
output_dir = ../runs/my-composite

# Matrix formulation; native defaults apply to omitted keys.
n1 = 128
m1 = 900
n2 = 32
functionality = 8
m4 = 6
density = 0.1
target_density = 0.8

# Particle geometry and interactions are configured independently.
np_count = 32
np_radius = 30
np_seed = 20261014
np_np_interaction = wca
np_matrix_interaction = wca
```

## Configure the particles

Lengths are angstrom, energies kcal/mol, and bead masses g/mol, using LAMMPS
`real` units. Numeric values must be finite and positive. NP keys reject
duplicates, including underscore/hyphen aliases.

| Key | Default | Meaning |
|---|---|---|
| `np_count` | Required | Integer particle count, 1..512. Use a plain native system for the zero-NP control. |
| `np_radius` | `30` | FCC carving radius measured to bead centers. |
| `np_representation` | `rigid_fcc` | Bead-built, translating and rotating rigid colloids; currently the supported representation. |
| `np_spacing` | `2^(1/6)*np_sigma` | FCC nearest-neighbor bead spacing. The radius must include at least the first shell. |
| `np_sigma` | `6.44584366` | NP bead contact parameter; used for NP–NP pairs and the default matrix mixing rule. |
| `np_epsilon` | `1.01287845` | NP bead reference well depth. |
| `np_bead_mass` | `74` | Each constituent bead's mass; independent of native matrix `mass`. |
| `np_seed` | `20261014` | Positive signed 32-bit seed for NP placement, orientations, matrix repacking and NP velocities. |
| `np_packing_min_distance` | `4.5` | Minimum intermolecular matrix bead separation during packing; covalent separations are preserved. |

The lattice constant is `sqrt(2)*np_spacing`. The local construction limit is
64 half-lattice indices along an axis and 2,000,000 total beads. Particles and
their contact clearances must fit the initial box, and NP pair cutoffs must be
less than half its length. Packing has bounded attempts and reports failure
instead of removing matrix molecules or accepting initial core overlaps.

At the defaults, radius 30 produces 429 beads per particle; the V22 example
with 122,430 matrix beads and 32 NPs has 136,158 beads and **10.0824042656 wt% NP**.
The default mapping is PDMS-equivalent mass and bead parameters. It is a rigid
colloid prototype rather than an atomistic silica force field or the integrated
LAMMPS `pair_style colloid` potential. Radius controls carving; sigma controls
constituent-bead contact distance.

## Configure interactions independently

| Key | Default | Meaning |
|---|---|---|
| `np_interaction` | `wca` | Both channels: `wca`, `lj`, or stage-matched `native_pdms`. |
| `np_np_interaction` | Inherits `np_interaction` | NP–NP constituent pairs. |
| `np_matrix_interaction` | Inherits `np_interaction` | NP pairs with every matrix type, including oil constituents. |
| `np_np_epsilon` | `np_epsilon` | Optional independent NP–NP well depth. |
| `np_matrix_epsilon` | Chemistry-specific mixing | Optional uniform NP–matrix well depth, overriding the mixing rule for all matrix types. |
| `np_matrix_sigma` | Chemistry-specific mixing | Optional uniform NP–matrix sigma. |
| `np_np_cutoff` | `2.5*sigma` for LJ | Explicit angstrom cutoff for NP–NP LJ. Rejected for WCA. |
| `np_matrix_cutoff` | `2.5*sigma` per pair for LJ | Explicit uniform angstrom cutoff for NP–matrix LJ. Rejected for WCA. |

Explicit channel choices override `np_interaction` regardless of config order.
WCA uses an energy-shifted LJ potential cut at `2^(1/6)*sigma`; it is purely
repulsive and reaches zero energy and force at the cutoff. LJ includes the
attractive tail and shifts energy to zero at its cutoff; an explicit LJ cutoff
must exceed the LJ minimum for every affected pair.

Without matrix overrides, the original 300 K PDMS/MPS reference cross-pairs
are retained at default NP epsilon/sigma. Changing `np_epsilon` scales each
reference cross-epsilon by `sqrt(np_epsilon / 1.01287845)`; changing `np_sigma`
adds half its difference from `6.44584366` to each cross-sigma. This preserves
the established DMS/MPS cross-interaction factors. NP–NP uses `np_sigma` and
`np_np_epsilon`. Every resolved pair is recorded in `composite.json`.

For `wca` and `lj`, NP parameters stay fixed through all matrix temperature stages. Native
matrix–matrix pairs, bonded terms and reaction rules are unchanged. Internal
nonbonded pairs within each rigid particle are excluded. For example, keeping
NP–NP repulsive while adding matrix adhesion requires:

```ini
np_np_interaction = wca
np_matrix_interaction = lj
np_matrix_epsilon = 0.8
```

Changing only interaction modes or epsilon leaves initial geometry unchanged.
Changing sigma can change spacing and core clearance; explicitly set
`np_spacing` if lattice spacing should remain fixed during a sigma study.

## Packing, integration and provenance

The engine generates the native matrix first and preserves its exact data,
input and metadata in `precursor/`. It unwraps molecules using their **bond
graphs**, preserving branches, rings, grafts and pendant beads. Whole molecules
are rotated and translated into a periodic box around independently oriented
NPs. Initial spatial placement changes; molecular identities, conformations,
functional-site counts and all covalent topology are preserved.

The initial box volume increases by the total mass ratio, maintaining the
native initial density. Matrix mass is summed from actual atom types, including
different MPS masses. `oil_wt` and `filler_wt` still describe the **matrix before
NP addition**; their percentages in the complete composite are correspondingly
lower. NP weight percentage is NP mass divided by total composite mass.

The next unused atom type identifies NP beads: type 6 for V22/V35 and type 4
for the generic elastomer. Each NP has its own molecule ID and is excluded
from crosslink creation. Rigid integration and ordinary matrix integration use
disjoint groups. Matrix NPT controls pressure using a whole-composite pressure
compute and a matrix-only temperature compute. Only flexible bead coordinates
are affinely rescaled; rigid centers are remapped by LAMMPS.

During native elastomer minimization, NP coordinates are held fixed. The native
minimizer and its tolerances are retained. Rigid state uses `reinit no`; after
each box-changing phase, one fixed-volume 5 fs step restores constituent
coordinates/images before stage output or rigid-fix redefinition. This adds
three steps for ordinary curing, four for conversion-controlled curing. The
update uses the **runtime timestep**, so an early conversion halt does not
invalidate the next run. Native curing limits and stage lengths are retained.

Every composite package contains:

| File | Purpose |
|---|---|
| `data.CASE`, `in.CASE` | Composite initial structure and complete preparation/cure/equilibration input. |
| `np-template.csv` | Canonical centered NP bead coordinates, independent of placed particle orientations. |
| `composite.json`, `nanoparticles.json`, `CASE.info` | Identical resolved composite metadata, counts, masses, pair table, centers and nested native matrix info. |
| `request.conf`, `generator.conf` | Requested composite config and exact native backend config. |
| `siliconelab.json` | Composite family, resolved base/backend, source provenance and file inventory. |
| `precursor/` | Native matrix data/input/info before augmentation. |

All original native configurations without NP settings continue through their
original generation paths. No-overwrite and failed-generation cleanup apply to
composites too. `--capabilities` provides machine-readable family/base/NP-model
discovery for agents; `--help composite` lists the NP keys and defaults.

## Local examples and validation

The reference recipes in `examples/composites/` include V22/WCA,
V35 with copolymer oil and attractive NP–matrix pairs, generic star elastomer,
conversion-controlled elastomer, and the original 32-NP V22 formulation through
the new selector. Their separate local output batch is `runs/composite-local/`.
See [the local workflow record](../workflows/composites.json).

CTest covers every base, all generic strand architectures, PDMS/PMPS/copolymer
oil, filler and conversion control, independent NP channels, custom bead
parameters/masses, radii/counts, topology/geometry/force-field preservation,
periodic clearance, reproducibility, invalid configs and existing-run safety.
The original native parity and analysis tests remain required. Local LAMMPS
preflight and short dynamics are recorded separately from production results.
Full cure convergence and calibrated mechanical response remain user tests.

The composite metadata format is `siliconelab-composite-info`. The generalized
structural analyzer supports composite validation and matrix-only basic/topology
analysis, including the legacy NP prototype format. See
[composite analysis scope and invocation](analyzer.md#composite-basic-and-topology-analysis).
Do not pair the preserved matrix metadata with full composite data. Composite
trajectory, Part 2/Part 3 and NP-obstacle primitive-path analysis remain unsupported.
The subset-temperature NPT warning is expected because pressure is explicitly
overridden with the whole-composite compute; image-flag warnings can arise from
new periodic crosslinks. Validation records identify these separately from NP
rigid-body image consistency.

## Native-PDMS interaction series and bulk tensile workflow

`np_interaction = native_pdms` copies the exact native matrix pair potential,
sigma and cutoffs in each hot/cold stage. `np_lambda` multiplies epsilon for
both NP–NP and NP–matrix pairs; matrix–matrix coefficients stay unchanged.
At lambda one, constituent interactions match native PDMS. NP beads remain
rigid and nonreactive. The mode rejects independent force-parameter overrides;
`np_native_stages` records the actual resolved stage coefficients. Placement
and initial coordinates are identical when only lambda changes.

`tensile_erate` (per fs) and `tensile_strain` must be supplied together and are
recorded in `siliconelab.json`. The bulk workflow reads these model.conf values
and generates independent x/y/z inputs using the exact cold formation force
field. The current series uses 2e-8/fs, strain one, 300 K, independent transverse
NPT at 1 atm, a 5 fs timestep and a final 5 fs fixed-volume consistency step.
Stress is tensile-positive: true stress is minus the full composite axial
pressure in MPa; engineering stress is true stress times instantaneous transverse
area divided by its initial value. No extra preloading relaxation or velocity
reset is added to the previously validated NP tensile protocol.

`scripts/bulk_workflow.py` saves the usual batch plan/state, model configs,
generated packages, formation and tensile stage folders. An explicit compute
preparation and `scripts/validate_bulk_workflow.py` preflight/smoke gate must pass
before explicit submission. Formation → x → y → z dependencies serialize each
case, while all three axes read the same undeformed final formation snapshot.
Seven cases therefore permit at most seven simultaneous batch jobs. Submission
checks expanded account queue elements against the eight-independent-job limit,
saves accepted IDs atomically and refuses duplicate or automatic MD resubmission.
Each stage pins the executor, checks inputs and parent contracts, and verifies
full duration/sampling, counts, hashes and rigid NP geometry on its compute node.
This simulation workflow does not enable unsupported composite Part 2/Part 3.
