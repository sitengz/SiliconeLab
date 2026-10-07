# Generator configuration reference

SiliconeLab prepares three classes of systems through one C++ executable:
`siliconelab_generator --config FILE`. This reference lists **every supported
configuration key** in the imported generators, including defaults, units,
valid choices, dependencies, derived quantities, and the existing upstream case
families. The small files in `examples/` are runnable preparation checks; they
are only a subset of this configuration space.

The reference describes the versions pinned in `upstream.json`. Backend defaults
are not the values in an upstream `model.conf`: a saved configuration can
explicitly override those defaults. An omitted optional key keeps the backend
default. “Unset” below means omit the key, rather than write an empty value.

## Contents

- [System choices and scientific scope](#system-choices-and-scientific-scope)
- [Common configuration and file syntax](#common-configuration-and-file-syntax)
- [Oil: complete settings](#oil-complete-settings)
- [Elastomer: complete settings](#elastomer-complete-settings)
- [Coating: complete settings](#coating-complete-settings)
- [Cross-system meanings and fixed simulation settings](#cross-system-meanings-and-fixed-simulation-settings)
- [Existing case families in the source repositories](#existing-case-families-in-the-source-repositories)
- [Generated files and provenance](#generated-files-and-provenance)
- [Source references and integration boundaries](#source-references-and-integration-boundaries)

## System choices and scientific scope

| Selection | Material and architecture | Composition choices | Geometry and generated workflow |
| --- | --- | --- | --- |
| `system = oil` | Uncrosslinked silicone chains; uniform length or an explicit length/count histogram. | Pure PDMS, pure PMPS, or PDMS/PMPS copolymer; random, alternating, or block sequence. | Initial periodic bulk, plus an input that converts the equilibrated bulk into a surface film. |
| `system = elastomer` | PDMS network precursors with linear, ring, star, or grafted strands and functional crosslinkers. | Functional-group stoichiometry; optional five-bead moderators and neutral PDMS filler; optional conversion-controlled curing. | Bulk or a directly generated wall-bounded film; main preparation input and final MSD sampling. |
| `system = coating`, `model = v22` | Straight, bifunctional PDMS strands; V22 formulation and reaction defaults. | No oil, PDMS, PMPS, or copolymer oil at a specified formulation weight fraction; optional moderators. | Bulk precursor or a directly generated film; films also get a dependent surface-measurement input. |
| `system = coating`, `model = v35` | Folded, bifunctional PDMS strands; V35 formulation and reaction defaults. | Same oil chemistry and loading controls as V22. | Same workflow choices as V22, with different strand geometry and placement capacity. |

Selecting `coating` alone does not select a film: omit `thickness` for a bulk
precursor, or supply a positive `thickness` for a film. Neither network generator
constructs the final reacted network during generation. Crosslinks form during
the generated LAMMPS calculation.

## Common configuration and file syntax

Start a job specification with its system, common output settings, and then the
scientific settings for that system:

```ini
system = coating
model = v22
output_dir = ../runs/my-coating

# Scientific details follow, using the tables below.
n1 = 128
m1 = 900
functionality = 8
oil = pdms
oil_length = 30
oil_wt = 10
```

### Common keys and commands

| Key or command | Type / choices | Default | Meaning |
| --- | --- | --- | --- |
| `system` | `oil`, `elastomer`, `coating` | Required | Selects the scientific backend. Values are lowercase and case-sensitive. |
| `model` | `v22`, `v35` | Required for coating | Selects the coating formulation. Rejected for oil and elastomer. |
| `output_dir` | Directory path | Required | New directory containing the complete run package. A relative path resolves beside the user's config file. |
| `--config FILE` | Command-line argument | Required for generation | Loads the common configuration; the file path resolves from the current working directory. |
| `--help` | Command-line argument | — | Shows the common interface. |
| `--help oil`, `--help elastomer`, `--help coating-v22`, `--help coating-v35` | Command-line arguments | — | Show the native scientific parameter help. |

The common executable accepts `--config FILE` or the help commands above. It
currently does **not** accept scientific command-line overrides such as
`--oil-wt 15`; edit the saved common config instead. Native backend help includes
options for the original standalone programs, whose command-line interfaces are
broader than the common entry point.

| Syntax rule | Behavior |
| --- | --- |
| Assignment | One `key = value` per line. Selection can technically appear anywhere; placing it first makes a config easier to read. |
| Whitespace and comments | Blank lines are ignored. `#` starts a comment, even inside quotes. |
| Quotes | Surrounding single or double quotes are stripped. Values are not shell expressions. |
| Key spelling | Underscores and hyphens are equivalent, e.g. `target_density` and `target-density`. A leading `--` on a key is accepted but unnecessary. |
| Repetition | Common keys cannot repeat. Native scalar settings generally take their last assigned value; explicit-setting conflicts still apply. Oil `chain_count` is intentionally repeatable. |
| Unknown keys | Rejected by the selected backend. Keys for one system are not automatically translated into another system's names. |
| Empty or nested input | Empty values, sections, YAML/JSON objects, and environment substitutions are unsupported. |
| Output naming | The last component of `output_dir` becomes the case name; use only ASCII letters, digits, `.`, `_`, and `-`. Parent directories may contain spaces. Output paths containing `#` or line breaks are unsupported. |
| Existing output | Every run directory must be new, including when an existing directory is empty. A failure removes only the new package created by that invocation; parent directories may remain. |
| Reserved native keys | `output` is replaced by `output_dir`. `config` and `help` cannot be config keys. |

Lengths and distances below are in **Å**, densities in **g/cm³**, and masses in
**g/mol** (numerically equivalent to bead masses in LAMMPS `real` units).
Percentages use the scale **0–100**, not 0–1. Counts and step intervals are
integers. Oil seeds require positive signed 32-bit integers. Elastomer and
coating seeds are parsed as signed 32-bit integers and then cast to unsigned
32-bit values; use positive values, especially for seeds also used by LAMMPS.
The tables describe supported scientific values rather than incidental parser
behavior for negative sentinel values.

## Oil: complete settings

Select `system = oil`. This backend has **14 scientific key spellings**,
including the `n` and `m` aliases. Its default chemistry is **pure PMPS**, not
PDMS: explicitly set `mps_percent = 0` for pure PDMS.

### Chain population and chemistry

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `length` | Positive integer | `16` | Repeat-unit positions per chain in the uniform-length mode. Not allowed with any `chain_count` row. |
| `n` | Positive integer | Same as `length` | Exact alias for `length`, including its conflicts. |
| `chains` | Positive integer | `625` | Number of chains in uniform-length mode. Not allowed with `chain_count`. |
| `m` | Positive integer | Same as `chains` | Exact alias for `chains`, including its conflicts. |
| `chain_count` | Two positive integers: `LENGTH COUNT`; repeatable | Unset | Explicit number of chains at each length. Lengths must be unique across rows. Cannot be combined with `length`, `chains`, `n`, or `m`. |
| `mps_percent` | Number, `0 <= X <= 100` | `100` | MPS fraction of repeat positions: `0` is PDMS, `100` is PMPS, and an intermediate value is copolymer. Mutually exclusive with `mps_wt`. |
| `mps_wt` | Number, `0 <= X <= 100` | Unset | MPS repeat-unit mass percentage within the oil, accounting for different DMS/MPS masses. Replaces `mps_percent`; do not explicitly set both. |
| `sequence` | `random`, `alternating`, `block` | `random` | Placement of the allocated MPS repeat positions within each chain. The chemistry allocation is performed before choosing their positions. |

### Geometry, packing, and randomness

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `density` | Positive number | `0.1` | Initial bulk mass density; determines the cubic box volume. |
| `target_density` | Positive number, at least `density` | `0.8` | Target density of the scripted 800 K compression. |
| `film_padding` | Positive number | Effective `20` | Requested extra space at each z face in the dependent film calculation. Values below 20 Å are raised to 20 Å; the runtime input may increase padding further to protect unwrapped chains. |
| `min_separation` | Number, `0 < X < 15` | `4.5` | Initial minimum intermolecular bead separation used during packing. |
| `seed` | Positive integer | `20260727` | Chain composition allocation, sequence construction, conformations, rotations, and packing randomness. |
| `velocity_seed` | Positive integer | `492845` | Seed for initial LAMMPS bulk velocities; independent of the structure-generation seed. |

### Oil combinations and derived quantities

| Scientific choice | Required settings or changes |
| --- | --- |
| Pure PDMS | `mps_percent = 0`; select `length`/`chains` or a `chain_count` histogram. |
| Pure PMPS | `mps_percent = 100`, or retain the backend chemistry default. |
| Copolymer specified by repeat fraction | Set `mps_percent` strictly between 0 and 100; choose a `sequence`. |
| Copolymer specified by mass fraction | Set `mps_wt`; omit explicit `mps_percent`; choose a `sequence`. |
| Uniform chain length | Set `length` and `chains`, or use their defaults/aliases. |
| Polydisperse population | Supply one `chain_count = LENGTH COUNT` row per distinct length. The input histogram, not a fitted distribution, defines the population. |

`random` chooses MPS positions randomly and reproducibly. `alternating` spreads
MPS positions as evenly as possible; it is not necessarily a strict 1:1 sequence
at other compositions. `block` places one central contiguous MPS block. For
pure PDMS or PMPS, the sequence setting does not change the chemistry.

A DMS repeat is one bead of mass 74.0. An MPS repeat consists of a backbone bead
of mass 59.1204 and a pendant bead of mass 77.106, totaling 136.2264. Thus a chain
with `N` repeats and `k` MPS repeats has:

```text
beads_per_chain = N + k
chain_mass = (N - k) * 74.0 + k * 136.2264
MPS_weight_percent = 100 * k * 136.2264 / chain_mass
```

The requested global MPS count is rounded to an integer. For uniform lengths,
neighboring per-chain MPS counts realize the global count. For unequal lengths,
the backend allocates counts proportionally using largest remainders. The native
`.info` records requested and realized composition, histograms, masses, and
atom/topology counts. System size is constrained by 32-bit LAMMPS atom IDs;
histogram chain/repeat totals are also checked against 32-bit limits.

There is no oil key named `oil`, `thickness`, `PDI`, or `mps_distribution`.
Oil geometry starts as bulk and its film input depends on the equilibrated
bulk snapshot; it does not construct a separate film from a thickness request.

## Elastomer: complete settings

Select `system = elastomer`. This backend has **27 scientific keys**.
Molecules are numbered in component order: strands, crosslinkers, moderators,
then filler. The defaults are a linear PDMS network precursor with no moderators
and no filler.

### Strand architecture and reactive groups

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `strand_topology` | `linear`, `ring`, `star`, `grafted` | `linear` | Strand architecture. Use `grafted`, not `graft`, `comb`, or `bottlebrush`, as the value. |
| `strand_length` | Positive integer | `128` for linear/ring; `32` per arm for star | Total strand beads for linear/ring; **beads per arm** for star. Explicitly forbidden for grafted strands. |
| `strand_count` | Positive integer | `900` | Number of strand molecules, regardless of architecture. |
| `strand_functionality` | Integer; architecture-dependent | `2` for linear/ring; derived for star/grafted | Number of reactive sites per strand. Linear strands must have 2. Ring strands allow 2 through their bead count. Star/grafted functionality is derived; an explicit value must match. |
| `strand_arm_count` | Integer: `3`, `4`, `6`, `8` | `4` when star | Number of star arms. Can be explicitly set only for `strand_topology = star`. |
| `backbone_length` | Positive integer | `64` when grafted | Grafted strand backbone beads. Explicit setting requires `grafted`. |
| `side_chain_length` | Positive integer | `16` when grafted | Beads in each grafted side chain. Explicit setting requires `grafted`. |
| `graft_spacing` | Nonnegative integer | `12` when grafted | Number of ungrafted backbone beads between grafts; the graft interval is `graft_spacing + 1`. Explicit setting requires `grafted`. |
| `graft_functional_fraction` | Number, `0 < X <= 100` | `40` when grafted | Requested percentage of side-chain ends that are reactive; rounded to whole ends, with at least one. Explicit setting requires `grafted`. |
| `strand_reactive_distribution` | `regular`, `random` | `regular` | Placement of ring reactive sites. `random` requires ring topology; other architectures use their defined endpoints/sites. |
| `strand_reactive_seed` | Integer seed | `20260810` | Random ring-site placement. Does not determine crosslinker reactive sites. |

### Crosslinkers, stoichiometry, and conversion

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `crosslinker_length` | Positive integer compatible with functionality | `32` | Beads per crosslinker chain. Its molecule count is derived, not supplied. |
| `functionality` | Integer, `3 <= F <= min(16, crosslinker_length)` | `4` | Reactive sites on each crosslinker; different from `strand_functionality`. |
| `stoichiometry` | `A:B`, with two positive integers | `1:1` | Ratio of strand functional groups to crosslinker functional groups, not a molecule-count ratio. Requires an exact whole-crosslinker count. |
| `crosslink_distribution` | `random`, `regular` | `random` | Initial placement of reactive sites along each crosslinker chain. |
| `crosslink_seed` | Integer seed | `20260722` | Random crosslinker-site selection. Not the runtime bond-creation RNG seed. |
| `target_conversion` | Number, `0 < X <= 100` | Unset | Enables conversion-controlled curing. Percentage of the stoichiometric maximum new bonds; the derived target must be at least one bond. Requires initial `density <= 0.5` and `target_density >= 0.5`. |

### Moderators and neutral PDMS filler

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `moderator_count` | Nonnegative integer | `0` | Optional five-bead star moderators: neutral center with four reactive arms. Their reactive sites are extra and excluded from the stoichiometric count and conversion target. Moderator length is fixed at five beads. |
| `filler_length` | Positive integer | Unset; no filler | Repeat/bead count per neutral PDMS filler chain. Must be supplied together with `filler_wt`. |
| `filler_wt` | Number, `0 < X < 100` | Unset; no filler | Filler percentage of complete-model mass. Filler chain count is derived and rounded. Must be supplied together with `filler_length`. |
| `filler_seed` | Integer seed | `20260727` | Filler conformations and packing randomness. Does not activate filler by itself. |
| `filler_min_separation` | Number, `0 < X < 15` when filler is used | `4.5` | Filler minimum bead distance from other components and other filler chains. |

### Mass, geometry, and initial state

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `mass` | Positive number | `74.0` | Common PDMS bead mass for strands, crosslinkers, moderators, and filler. |
| `density` | Positive number | `0.1` | Initial bulk density, or film density based on nominal wall-free material volume. |
| `target_density` | Positive number | `0.8` | Scripted compression target using the same density convention. Conversion mode adds the 0.5 g/cm³ constraint above. |
| `thickness` | Positive number to request film | Unset; bulk | **Nominal 300 K wall-free material thickness**, not the full simulation-box height. Zero is rejected; omitted/negative native sentinel selects bulk. |
| `seed` | Integer seed | `5489` | Native placement/moderator randomness and initial LAMMPS velocity seed. |

### Architecture rules and derived sizes

| Architecture | Length, reactive-site, and compatibility rules |
| --- | --- |
| Linear | `strand_length` is the full chain length; both ends are reactive; `strand_functionality = 2`. Do not supply star or graft controls. |
| Ring | At least four beads; cyclic bonds/angles/dihedrals. `2 <= strand_functionality <= strand_length`. Regular reactive placement requires `strand_length` divisible by `strand_functionality`; random placement samples distinct sites. |
| Star | `strand_length` is the arm length, default 32; all outer arm ends are reactive. Functionality equals arm count. No graft controls. |
| Grafted, comb-like | `strand_topology = grafted` with positive `graft_spacing`; the backbone ends are neutral and selected side-chain ends are reactive. |
| Grafted, dense bottlebrush | `strand_topology = grafted`, `graft_spacing = 0`; one side chain on each backbone bead. This is a spacing choice, not a separate topology keyword. |

Star total beads and functionality are derived as follows:

| `strand_arm_count` | Neutral center beads | Total beads per strand | Derived `strand_functionality` |
| --- | --- | --- | --- |
| `3` | `1` | `3 * strand_length + 1` | `3` |
| `4` | `1` | `4 * strand_length + 1` | `4` |
| `6` | `2` | `6 * strand_length + 2` | `6` |
| `8` | `3` | `8 * strand_length + 3` | `8` |

For a grafted strand:

```text
side_chain_count = ceil(backbone_length / (graft_spacing + 1))
strand_beads = backbone_length + side_chain_count * side_chain_length
strand_functionality = clamp(round(side_chain_count * graft_functional_fraction / 100),
                             1, side_chain_count)
```

Functional graft ends are selected regularly along the graft sequence; random
ring-site selection does not apply to grafted strands. The requested and
realized functional fractions are recorded. Derived star/grafted strand bead
counts above 100,000,000 are rejected by the backend.

### Exact functional-group stoichiometry

For `stoichiometry = A:B`, with resolved strand functionality `Fs` and
crosslinker functionality `Fc`:

```text
strand_groups = Fs * strand_count
crosslinker_count = strand_groups * B / (Fc * A)
crosslinker_groups = Fc * crosslinker_count
maximum_new_bonds = min(strand_groups, crosslinker_groups)
```

`crosslinker_count` must be an exact integer and no larger than 100,000,000.
For default linear strands and `A:B = 1:1`, it is `strand_count / 2`:
900 strands imply 450 crosslinkers, rather than 900. Moderator sites are
excluded. Regular crosslinker sites occupy bead positions 1, 3, 5, ... and
require `functionality <= ceil(crosslinker_length / 2)`. Random placement chooses
distinct sites using `crosslink_seed`.

For filler loading `w = filler_wt / 100`:

```text
base_mass = mass of strands + crosslinkers + moderators
requested_filler_mass = base_mass * w / (1 - w)
filler_chain_mass = filler_length * mass
filler_count = max(1, round(requested_filler_mass / filler_chain_mass))
```

Filler remains neutral and is packed in the central 40% of the initial z box,
from `-0.20 * Lz` to `+0.20 * Lz`. The realized weight fraction can differ from
the request because molecule counts are integers.

### Conversion-controlled and film modes

With an explicit `target_conversion = X`:

```text
target_new_bonds = floor(maximum_new_bonds * X / 100)
```

The generated input equilibrates at 800 K, compresses without reactions to
0.5 g/cm³, and cures there with bond-creation probability 0.5 for at most
5,000,000 steps. It checks the bond target each timestep, then removes reaction
fixes, compresses to `target_density`, equilibrates, cools, and performs final
sampling. A plateau below target is possible; simultaneous final-step reactions
can overshoot the integer target. The config records an intended target, not a
measured achieved conversion. Omitting the key keeps the native fixed-duration
curing schedule with probability 0.1.

For a film, the backend constructs:

```text
material_volume = total_mass / (density * 0.602)
Lz = thickness + 2 * cold_wall_cutoff
Lx = Ly = sqrt(material_volume / thickness)
```

`cold_wall_cutoff` is the model's 300 K repulsive-wall cutoff. Compression is
lateral at fixed box `Lz`; material-density settings use the nominal wall-free
volume. The `.info` records nominal thickness, total box height, wall cutoff,
and density convention. This meaning differs from the coating `thickness` key.
There are no elastomer keys for `oil`, PMPS/copolymer chemistry,
`bond_length`, `spacing`, or the coating `surface_*` controls.

## Coating: complete settings

Select `system = coating` and `model = v22` or `v35`. The two formulations
share **28 scientific keys**. The base network is PDMS; PDMS/PMPS chemistry
controls apply to the oil component, not to the network strands.

### Formulation defaults

| Default or fixed model choice | V22 | V35 |
| --- | --- | --- |
| `n1`: strand beads | `128` | `384` |
| `m1`: strand molecules | `900` | `306` |
| Strand shape | Straight linear | Folded serpentine |
| Strand functionality | Fixed `2`, reactive ends | Fixed `2`, reactive ends |
| `n2`: crosslinker beads | `32` | `32` |
| `functionality`: crosslinker sites | `8` | `4` |
| Derived `m2` at these defaults | `225` | `153` |
| Oil component | Absent | Absent |
| `n4`, `m4`: moderators | `5`, `6` | `5`, `6` |
| Geometry | Bulk | Bulk |
| Fixed runtime bond-creation probability | `0.1` | `0.5` |

### Network, crosslinkers, and moderators

| Key | Type / valid values | Default: V22 / V35 | Meaning and dependencies |
| --- | --- | --- | --- |
| `n1` | Positive integer when `m1 > 0`; nonnegative otherwise | `128` / `384` | Beads in each bifunctional PDMS network strand. |
| `m1` | Nonnegative integer | `900` / `306` | Number of network strands. Exact crosslinker stoichiometry must be possible. |
| `n2` | Positive integer compatible with functionality | `32` / `32` | Beads per crosslinker chain. |
| `functionality` | Integer, `3 <= F <= min(16, n2)` | `8` / `4` | Reactive sites per crosslinker; also determines its derived molecule count. |
| `crosslink_distribution` | `random`, `regular` | `random` / `random` | Placement of reactive sites along crosslinker chains. Regular sites are 1, 3, 5, ... and require `F <= ceil(n2 / 2)`. |
| `crosslink_seed` | Integer seed | `20260722` / `20260722` | Random crosslinker reactive-site selection. |
| `n4` | Nonnegative integer; must be `5` when `m4 > 0` | `5` / `5` | Moderator bead count. The implemented moderator is a neutral center plus four reactive arms. |
| `m4` | Nonnegative integer | `6` / `6` | Moderator molecule count; `0` disables moderators. Unlike elastomer, moderators are on by default. |

### Oil chemistry and loading

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `oil` | `pdms`, `pmps`, `copolymer` | Unset: no oil | Explicitly enables the selected oil component. Literal `oil = none` is rejected; omit the key for the no-oil control. |
| `oil_length` | Positive integer; at least 2 for copolymer | Unset | Repeat-unit positions per oil chain. Required with `oil`, together with `oil_wt`. |
| `oil_wt` | Number, `0 < X < 100` | Unset | Oil percentage of total formulation mass, including network, crosslinkers, moderators, and oil. Required with `oil`. |
| `mps_percent` | Number, `0 < X < 100` | Unset | MPS repeat-position percentage within copolymer oil. Required unless using `mps_wt`; rejected for pure PDMS/PMPS. |
| `mps_wt` | Number, `0 < X < 100` | Unset | MPS repeat-unit mass percentage within copolymer oil. Alternative to `mps_percent`; not the oil loading in the full coating. |
| `mps_distribution` | `fixed`, `balanced` | `fixed` | Copolymer MPS-count allocation across oil chains. `balanced` requires `oil = copolymer`. |
| `sequence` | `random`, `alternating`, `block` | `random` | MPS positions within copolymer chains. Explicit sequence controls are rejected for pure PDMS/PMPS. |
| `oil_seed` | Integer seed | `20260727` | Oil sequence, conformations, and packing randomness. Setting it alone does not enable oil. |
| `oil_min_separation` | Number, `0 < X < 15` when oil is enabled | `4.5` | Minimum oil/moderator-to-previous-component bead distance during initial placement; also used for moderator overlap rejection. |

### Mass, placement, and geometry

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `mass` | Positive number | `74.0` | Common DMS bead mass. MPS backbone and pendant masses remain 59.1204 and 77.106. |
| `density` | Positive number | `0.1` | Initial total mass density; determines initial box volume. |
| `target_density` | Positive number | `0.8` | Scripted compression target. |
| `bond_length` | Positive number | `2.801` | Initial network/crosslinker/moderator construction bond length. Does not rewrite bonded force-field coefficients or the shared oil-component geometry. |
| `spacing` | Positive number | `7.0` | Initial formulation strand/crosslinker placement spacing; affects folded-network capacity. |
| `thickness` | Positive number to request film | Unset; bulk | Initial film simulation-box height `Lz`, fixed during the main film curing input. Zero is rejected; omitted/negative native sentinel selects bulk. |
| `seed` | Integer seed | `5489` | Star-moderator placement seed; not the oil seed or runtime bond-creation seed. |

### Dependent surface stage

| Key | Type / valid values | Default | Meaning and dependencies |
| --- | --- | --- | --- |
| `surface_padding` | Finite number, greater than `15` | `50` | Added vacuum at each z face for the dependent surface-measurement stage. |
| `surface_relax_steps` | Positive integer | `10000000` | 300 K surface NVT relaxation steps. At fixed 5 fs timestep, the default is 50 ns. |
| `surface_production_steps` | Positive integer | `10000000` | 300 K surface-pressure production steps; default 50 ns. |
| `surface_sample_every` | Positive integer | `1000` | Pressure sampling interval in MD steps; default 5 ps. |

Surface settings are validated even for bulk configs, but the separate surface
input/scripts are produced only when `thickness > 0`. Sampling intervals are not
automatically adjusted to give a desired number of samples: choose them relative
to both step counts.

### Valid oil combinations and rounding

| Oil case | Settings to use | Settings to omit |
| --- | --- | --- |
| No-oil control | Select `model`, network, geometry, and moderator controls. | Omit `oil`, `oil_length`, `oil_wt`, `mps_percent`, `mps_wt`, and explicit `sequence`. Oil seeds/separation/distribution alone do not enable oil. |
| PDMS oil | `oil = pdms`, positive `oil_length`, `0 < oil_wt < 100`. | MPS composition and explicit sequence keys; do not set `oil_wt = 0` for a no-oil control. |
| PMPS oil | `oil = pmps`, positive `oil_length`, `0 < oil_wt < 100`. | MPS composition and explicit sequence keys. |
| Copolymer by repeat percentage | `oil = copolymer`, `oil_length >= 2`, `oil_wt`, `mps_percent`, optionally sequence/distribution. | `mps_wt`. |
| Copolymer by repeat mass percentage | `oil = copolymer`, `oil_length >= 2`, `oil_wt`, `mps_wt`, optionally sequence/distribution. | `mps_percent`. |

Pure PDMS/PMPS oils use the endpoints 0/100 internally. Unlike the standalone
oil generator, coating copolymer composition controls require a strict interior
percentage; use `oil = pdms` or `oil = pmps` for pure oils.

For every coating, `m2 = 2 * m1 / functionality` must be an exact integer and
no larger than 100,000,000. There is no configurable `stoichiometry` or
`strand_functionality` for coatings. Moderator reactive sites are additional
and excluded from this count. The complete formulation must contain beads;
zero component counts are accepted only when the remaining composition and
placement are valid.

With `w = oil_wt / 100`:

```text
base_mass = mass * (n1 * m1 + n2 * m2 + n4 * m4)
requested_oil_mass = base_mass * w / (1 - w)
oil_chain_mass = (oil_length - k) * mass + k * 136.2264
m3 = max(1, round(requested_oil_mass / estimated_oil_chain_mass))
```

`k` denotes MPS repeats per chain. With `mps_distribution = fixed`, the requested
per-chain MPS count is rounded and clamped to at least one of each repeat type
for copolymers, then every oil chain uses that count. With `balanced`, the
backend estimates chain count from the unrounded requested composition, rounds
the **global** MPS count, and spreads it across neighboring per-chain counts.
At extreme fractions and short lengths, individual balanced chains may become
pure even though the requested population is a copolymer. Inspect `.info` for
realized oil weight fraction and MPS repeat composition.

Oil is initially confined to `-0.20 * Lz` through `+0.20 * Lz`, the central 40%
of the box. Moderators occupy an upper-middle band, `+0.28 * Lz` through
`+0.38 * Lz`, and avoid previously placed components. These bands are fixed,
not configurable concentration profiles. V22 constructs straight strands from
the bottom upward; crosslinkers are placed from the top downward. V35 folds
its network strands. V35 requires adequate lateral space for the folded
footprint and enough row/layer capacity for `m1`; a nominally valid parameter
set may still fail geometric placement or overlap checks.

### Coating bulk and film construction

```text
volume = total_mass / (density * 0.602)
bulk: Lx = Ly = Lz = cbrt(volume)
film: Lz = thickness; Lx = Ly = sqrt(volume / thickness)
```

The main film input uses separate repulsive z guard walls, nonperiodic z, and
lateral compression/pressure control. Its surface input reads the **film's**
final `data.<case>.npt_eq`, adds `surface_padding` at each z face, and records
surface-pressure production after relaxation. It does not turn an already
reacted periodic bulk network into the independently cured film counterpart.
The upstream film-series helper measures bulk `Lz` to choose a new film's
thickness; that helper is not yet a SiliconeLab command.

### Rejected coating inputs from older conventions

| Rejected key | Reason / replacement |
| --- | --- |
| `m2` | Derived exactly from `2 * m1 / functionality`. |
| `n3`, `m3` | Legacy direct oil-size/count controls; use `oil`, `oil_length`, and `oil_wt`. |
| `filler_length`, `filler_wt` | Obsolete coating controls; use the explicit oil controls. These names remain valid for elastomer's neutral PDMS filler. |
| `oil = none` | Omit `oil` and chemistry/loading controls to request no oil. |
| `target_conversion` | Implemented only by the elastomer backend, not by V22/V35. |

## Cross-system meanings and fixed simulation settings

### Similar names with different meanings

| Concept | Oil | Elastomer | Coating |
| --- | --- | --- | --- |
| Chain/strand length | `length`: repeat positions; an MPS repeat adds a pendant bead. | `strand_length`: full linear/ring size or per-arm star size; graft size is derived. | `n1`: PDMS strand beads; `oil_length`: oil repeat positions. |
| Oil loading | The entire system is oil; `mps_wt` sets chemistry within it. | `filler_wt`: neutral PDMS filler fraction of complete-model mass. | `oil_wt`: oil fraction of complete-coating mass; `mps_wt`: composition within that oil. |
| Crosslinker functionality | No network crosslinkers. | `functionality`, distinct from `strand_functionality`. | `functionality`; strand functionality is fixed at 2. |
| Film thickness | No thickness input; film comes from the equilibrated bulk. | `thickness`: nominal wall-free material thickness; the box includes wall offsets. | `thickness`: initial film box `Lz`; later surface stage adds padding. |
| Geometry controls | `min_separation`, `film_padding`. | `filler_min_separation`; backbone bond length and spacing are fixed. | `bond_length`, `spacing`, `oil_min_separation`, `surface_padding`. |
| Default optional components | Oil only; PMPS chemistry by default. | Moderators/filler off. | Oil off; six moderators on. |
| Conversion target | Not applicable. | `target_conversion`, in percent. | No target key; model-specific fixed-duration curing. |
| Seed roles | `seed` controls structures; `velocity_seed` controls velocities. | `seed` also controls initial velocities; site/filler seeds are separate. | `seed` is for moderators; `oil_seed` and `crosslink_seed` are separate. |

### Settings compiled into the current backends

The keys listed above are the configuration interface. The following scientific
and execution settings are **not additional accepted config keys**:

| Setting | Current behavior |
| --- | --- |
| Temperature / timestep / force field | Generated calculations use the upstream 800 K preparation and 300 K final stages, 5 fs timestep, and explicit model coefficients. There is no common `temperature`, `timestep`, `force_field`, or arbitrary coefficient override. |
| Oil masses and initial bonded geometry | DMS mass 74.0; MPS backbone 59.1204 and pendant 77.106. The standalone oil has no `mass` or `bond_length` key. |
| Elastomer construction geometry | PDMS bond length 2.801 Å and initial placement spacing 7.5 Å are fixed. |
| DMS/MPS mixing in oil/coating | Cross chemistry uses `0.579966 * geometric epsilon` and arithmetic sigma; model pair matrices and bonded coefficients are fixed. `mass` in coating changes DMS mass, not these coefficients. |
| Nonbonded 1–4 convention | The current generators write `special_bonds lj 0 0 0.5`. |
| Runtime bond-creation random seed | Network/coating `fix bond/create` uses fixed seed `348154`; `crosslink_seed` instead chooses the initial crosslinker sites. |
| Runtime bond-creation probability | Elastomer fixed-duration mode: 0.1; elastomer target mode: 0.5; V22: 0.1; V35: 0.5. There is no config key to set the probability independently. |
| Film wall conventions | Native z-wall models/cutoffs are fixed; oil film padding and coating surface padding are the exposed controls. |
| Cluster resources | Inherited Slurm templates contain Nova partition/module/resource assumptions. Host, account, queue, walltime, node/task count, memory, and mail settings are not generator config keys. Cluster profiles are a later integration step. |
| Analysis and requested properties | Analyzer windows, block sizes, convergence criteria, properties, and tensile/layer-dynamics jobs are outside this generator config. They belong to later simulation/analysis specifications. |

### Generated simulation budgets

These are planned input schedules; local generation does not execute them or
establish equilibration, target conversion, or convergence.

| Backend / stage | Default or fixed budget |
| --- | --- |
| Oil bulk | 11,000,000 steps: 800 K relaxation/compression stages, cooling, then 5,000,000 final 300 K NPT steps. No crosslinking. |
| Oil film, after bulk | 100,000 initial guarded lateral-relaxation steps, 10,000,000 NVT relaxation steps, 10,000,000 NVT pressure-production steps. These counts are fixed; `film_padding` is configurable. |
| Elastomer fixed-duration preparation | 7,000,000 main steps, followed by 1,000,000 MSD-production steps. |
| Elastomer conversion-controlled preparation | At most 11,000,000 main steps, including a curing hold of at most 5,000,000 steps; followed by 1,000,000 MSD-production steps. Early target detection can shorten curing. |
| V22/V35 main preparation | 7,000,000 main steps, followed by 1,000,000 MSD-production steps; the main stage durations are fixed. |
| V22/V35 film surface stage | `surface_relax_steps` plus `surface_production_steps`; defaults 10,000,000 + 10,000,000, sampled every `surface_sample_every` steps. |

The independent network/coating MSD production writes timestep-zero plus every
1,000 steps through 1,000,000: 1,001 expected frames. Those trajectory settings
are fixed. At 5 fs, 1,000,000 steps is 5 ns. Full physical model definitions and
stage-specific LAMMPS instructions remain in the pinned sources linked below.

## Existing case families in the source repositories

This catalog gives the broader experimental/design space already represented
in the three repositories. These source cases are **not all imported as runnable
SiliconeLab configs**. To port a case, preserve its scientific settings, add
`system` and coating `model` when needed, and replace its native `output` with a
new `output_dir`. The small examples in SiliconeLab are independent size-reduced
preparation checks, not replacements for these production cases.

### Seven oil chemistries and the corresponding coating series

The oil repository has seven numbered configurations. The coating repository
uses the same oil chemistries at 5 and 10 wt% under both V22 and V35.

| Upstream oil case | Chemistry / sequence | Oil repeat length | MPS repeat percentage | Chains in the oil-only case | Coating choices |
| --- | --- | --- | --- | --- | --- |
| `01` | PDMS | `30` | `0` | `3333` | `oil = pdms`, `oil_length = 30`, `oil_wt = 5` or `10`. |
| `02` | PMPS | `12` | `100` | `8333` | `oil = pmps`, `oil_length = 12`, `oil_wt = 5` or `10`. |
| `03` | Random copolymer | `179` | `5` | `559` | `oil = copolymer`, length 179, `mps_percent = 5`, loading 5 or 10 wt%. |
| `04` | Random copolymer | `42` | `10` | `2381` | Copolymer, length 42, 10% MPS, loading 5 or 10 wt%. |
| `05` | Random copolymer | `44` | `10` | `2273` | Copolymer, length 44, 10% MPS, loading 5 or 10 wt%. |
| `06` | Random copolymer | `65` | `10` | `1538` | Copolymer, length 65, 10% MPS, loading 5 or 10 wt%. |
| `07` | Random copolymer | `15` | `50` | `6667` | Copolymer, length 15, 50% MPS, loading 5 or 10 wt%. |

The numbered oil configs target about 100,000 repeat positions. Their chain
counts **must not be copied into coatings**: coating `m3` is derived from the
complete formulation's mass and oil loading. The coating series uses
`mps_distribution = balanced` for copolymers, `sequence = random`, the V22 base
`n1 = 128`, `m1 = 900`, `functionality = 8`, or the V35 base `n1 = 384`,
`m1 = 306`, `functionality = 4`. Each formulation has 14 oil-containing cases
plus a `NoOil_0wt` control, totaling 30 bulk configurations. The upstream film
helper prepares corresponding independent films after completed bulk runs.

Sources: [oil series](https://github.com/sitengz/Silicone_Oil/blob/66b00f69b09736e607a18303c4619d5a6c089dc5/simulations/README.md),
[coating series](https://github.com/sitengz/Silicone_Coating/blob/169ad4f8a3cb5bb0d205c576e1a1897bd2d63430/simulations/README.md).

### Oil PDMS chain-length / PDI series

| Family | Cases | Actual generator input |
| --- | --- | --- |
| Monodisperse PDMS | Lengths `4`, `8`, `16`, `32`, `64`, `128`, all PDI 1; six cases. | Explicit `chain_count` histogram with a single occupied length; `mps_percent = 0`. |
| Polydisperse PDMS | Number-average repeat lengths 16, 32, 64; target PDIs 1.05, 1.10, 1.15, 1.20, 1.25, 1.30; eighteen cases. | Explicit integer `chain_count` rows constructed by the upstream series-building script. |

Together these are 24 `N*_PDI*` cases. They have at most 100,000 DMS beads and
exact integer-repeat number-average length. The N16/PDI1.30 histogram uses
lengths 4–34; other N16 cases are capped at 32, N32 cases at 64, and N64 cases
at 128. The generator has no `PDI`, `Mn`, `Mw`, or distribution-fitting key.
For pure PDMS histograms, using counts `c_i` at lengths `N_i`:

```text
number_average_repeat_length = sum(c_i * N_i) / sum(c_i)
weight_average_repeat_length = sum(c_i * N_i^2) / sum(c_i * N_i)
PDI = weight_average_repeat_length / number_average_repeat_length
```

Multiply these length averages by 74.0 to obtain the model's chain molar-mass
averages. Upstream `build_pdms_series.py` prepares the histograms and records
their realized properties; the unified generator consumes the resulting rows.
See the [PDMS case catalog](https://github.com/sitengz/Silicone_Oil/blob/66b00f69b09736e607a18303c4619d5a6c089dc5/simulations/README.md)
and [realized histogram tables](https://github.com/sitengz/Silicone_Oil/blob/66b00f69b09736e607a18303c4619d5a6c089dc5/simulations/pdms_series_tables.md).

### Elastomer architecture examples

The upstream ten-example collection covers more than the linear example in
SiliconeLab. All use 32-bead, four-functional crosslinkers at `stoichiometry = 1:1`
and disable moderators/filler unless explicitly edited.

| Upstream example | Strand-specific settings | Derived beads / functionality per strand |
| --- | --- | --- |
| `01_default` | Linear, `strand_length = 128`, `strand_count = 900`. | `128` beads, `2` sites. |
| `02_ring_bifunctional` | Ring, length 128, count 900, functionality 2, regular sites. | `128` beads, `2` sites. |
| `03_ring_tetrafunctional` | Ring, length 128, count 900, functionality 4, random sites, seed 20260810. | `128` beads, `4` sites. |
| `04_star_3arm` | Star, arm length 32, count 1200, 3 arms. | `97` beads, `3` sites. |
| `05_star_4arm` | Star, arm length 32, count 900, 4 arms. | `129` beads, `4` sites. |
| `06_star_6arm` | Star, arm length 32, count 600, 6 arms. | `194` beads, `6` sites. |
| `07_star_8arm` | Star, arm length 32, count 450, 8 arms. | `259` beads, `8` sites. |
| `08_grafted_comb` | Grafted, backbone 64, side chain 8, spacing 12, functional fraction 40%, count 1200. | `104` beads, `2` sites. |
| `09_grafted_bottlebrush` | Grafted, backbone 24, side chain 5, spacing 0, functional fraction 25%, count 750. | `144` beads, `6` sites. |
| `10_default_film` | Linear counterpart of example 01; the upstream run procedure supplies a measured film thickness. | `128` beads, `2` sites; an actual positive `thickness` is needed to select film. |

Additional knobs allow off-stoichiometric functional-group ratios, regular or
random crosslinker sites, moderators, filler loading, and conversion targets.
They are not distinct `system` values. See the [architecture discussion](https://github.com/sitengz/PDMS_Elastomer/blob/b67c4cea77893cc9eae41a671193c8b6c201606d/README.md)
and [original example configurations](https://github.com/sitengz/PDMS_Elastomer/tree/b67c4cea77893cc9eae41a671193c8b6c201606d/examples).

## Generated files and provenance

For `output_dir = ../runs/my-case`, all systems share this package convention:

| File | Purpose |
| --- | --- |
| `data.my-case` | Initial LAMMPS structure and topology. |
| `in.my-case` | Native bulk/network preparation input. |
| `submit.my-case.sh` | Inherited cluster submission template. |
| `my-case.info` | Native resolved scientific metadata: defaults, composition, architecture, geometry, seeds, counts, and planned simulation settings. |
| `request.conf` | Exact copy of the common configuration supplied by the user. |
| `generator.conf` | Scientific settings passed to the selected backend with an absolute native output path. |
| `siliconelab.json` | Schema version, `status = generated`, system/model, upstream repository/commit, data/info filenames, and package inventory. |

| System / condition | Additional generated files and dependencies |
| --- | --- |
| Oil, every config | `in.my-case.film`, `submit.my-case.film.sh`, `submit.my-case.pair.sh`. The dependent film input requires `data.my-case.npt_eq` from a successful bulk calculation. |
| Elastomer | Main input includes final MSD production. There is no separate coating-style surface input in this backend. |
| Coating, bulk | Main input includes final MSD production; no surface-stage package. |
| Coating, film | `in.my-case.surface`, `submit.my-case.surface.sh`, `submit.my-case.chain.sh`. Surface measurement requires the film's own completed `data.my-case.npt_eq`; the inherited chain launcher encodes the Slurm dependency. |

Manifest inventory paths are relative to the package. `generator.conf` records
supplied scientific settings, not an expansion of every default. Consult the
native `.info` for resolved and realized values. Completion of generation means
inputs were written, not that a simulation or property calculation succeeded.
Moving a package keeps its relative data/input references, but its saved
`generator.conf` still contains the original absolute generation output path;
edit the common request's `output_dir` to regenerate elsewhere.

## Source references and integration boundaries

The parameter inventory above is checked against the actual `apply_option`
handlers, settings defaults, architecture/composition resolvers, validation,
and output writers in these imported sources:

| Backend | Imported source | Source revision |
| --- | --- | --- |
| Oil | [oil generator](../vendor/Silicone_Oil/Generator/oil_generator.cpp) | `66b00f69b09736e607a18303c4619d5a6c089dc5` |
| Elastomer | [elastomer generator](../vendor/PDMS_Elastomer/Generator/pdms_elastomer_generator.cpp) | `b67c4cea77893cc9eae41a671193c8b6c201606d` |
| Coating V22 | [V22 formulation wrapper](../vendor/Silicone_Coating/V22/v22_generator.cpp) | `169ad4f8a3cb5bb0d205c576e1a1897bd2d63430` |
| Coating V35 and shared implementation | [V35 generator](../vendor/Silicone_Coating/V35/v35_generator.cpp), [shared oil component](../vendor/Silicone_Coating/Formulation/silicone_oil_component.hpp) | `169ad4f8a3cb5bb0d205c576e1a1897bd2d63430` |

For additional scientific context, see the pinned upstream
[oil README](https://github.com/sitengz/Silicone_Oil/blob/66b00f69b09736e607a18303c4619d5a6c089dc5/README.md),
[elastomer README](https://github.com/sitengz/PDMS_Elastomer/blob/b67c4cea77893cc9eae41a671193c8b6c201606d/README.md),
[V22 README](https://github.com/sitengz/Silicone_Coating/blob/169ad4f8a3cb5bb0d205c576e1a1897bd2d63430/V22/README.md),
[V35 README](https://github.com/sitengz/Silicone_Coating/blob/169ad4f8a3cb5bb0d205c576e1a1897bd2d63430/V35/README.md),
and [shared formulation description](https://github.com/sitengz/Silicone_Coating/blob/169ad4f8a3cb5bb0d205c576e1a1897bd2d63430/Formulation/README.md).

Source files are imported without edits; SHA-256 values and original paths are
recorded in `upstream.json`, and the original MIT licenses are preserved.
CMake compiles each backend separately, changes its entry-point name, and links
it into the common executable. No network access is needed to build.

Other tools in the upstream repositories are not yet part of this entry point:
the oil PDI-series constructor, the elastomer `tensile_test_generator` and
`layer_dynamics_generator`, the coating `generate_films.py` helper, and all
analyzers. Their settings should not be added to a generator config expecting
the common executable to interpret them.

Local parity checks cover the included small examples and several advanced
combinations. They do not prove convergence or validate every configuration in
this reference. Geometric placement, realizable composition, integer counts,
and native validation still constrain combinations within the listed ranges.
Cluster connection, job submission, MD execution, retrieval, and analysis remain
subsequent workflow steps. Review inherited cluster scripts for the intended
machine before submission.
