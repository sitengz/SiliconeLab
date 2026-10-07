# Common generator configuration

Run `siliconelab_generator --config FILE`. The configuration is a flat
`key = value` file, using the format of the source repositories. Put the system
selection first, then common output settings, followed by scientific details.

| Common key | Meaning |
| --- | --- |
| `system` | Required: `oil`, `elastomer`, or `coating` (lowercase). |
| `model` | Required for coating: `v22` or `v35`. Not accepted for other systems. |
| `output_dir` | Required: a new directory for the entire run package. Relative paths resolve beside the configuration file. |

Common settings cannot repeat. Scientific settings are passed to the selected
backend in their original order, including repeated oil `chain_count` rows.
Underscores and hyphens in keys are equivalent. Values can have surrounding
single or double quotes. `#` starts a comment, including inside quoted values,
as in the native config format. Output paths containing `#` or line breaks are
unsupported. There are no sections, environment substitutions, or YAML parsing.

The last component of `output_dir` becomes the case name and must contain only
ASCII letters, digits, dots, underscores, and hyphens. Parent directories can
contain spaces. This keeps the inherited LAMMPS filename tokens valid.

Native `output` is replaced by `output_dir`; do not specify both. Existing run
directories are rejected, including empty directories. On a generation error,
the new package is removed; its diagnostic is written to stderr. Parent
directories may remain. No existing package is overwritten or removed.

## Oil

```ini
system = oil
output_dir = ../runs/oil-N16

length = 16
chains = 12
mps_percent = 0
sequence = random
density = 0.1
target_density = 0.8
seed = 20260727
velocity_seed = 492845
```

`mps_percent = 0` selects PDMS repeats; `100` selects PMPS repeats. Intermediate
values select copolymer composition. The upstream default is **100**, so specify
the composition deliberately. `mps_wt` is an alternative weight-percent control
and cannot be combined with `mps_percent`.

For polydispersity, replace `length` and `chains` with repeated rows:

```ini
chain_count = 8 4
chain_count = 16 4
```

Each row is repeat length followed by molecule count. Do not combine the rows
with `length` or `chains`. Other controls include `sequence` (random, alternating,
block), `min_separation`, and `film_padding`. The package includes a bulk input
and a dependent film input. The film stage needs the equilibrated data produced
by the bulk LAMMPS calculation.

## Elastomer

```ini
system = elastomer
output_dir = ../runs/elastomer-linear

strand_topology = linear
strand_length = 16
strand_count = 12
strand_functionality = 2
crosslinker_length = 16
functionality = 4
stoichiometry = 1:1
moderator_count = 0
crosslink_distribution = random
crosslink_seed = 20260722
density = 0.1
seed = 5489
```

`strand_topology` supports linear, ring, star, and graft architectures. Their
architecture-specific constraints are still enforced by the original generator.
Controls include `strand_arm_count`, `backbone_length`, `side_chain_length`,
`graft_spacing`, `graft_functional_fraction`, `strand_reactive_distribution`,
and `strand_reactive_seed`.

Use `filler_length` with `filler_wt` for neutral PDMS filler; its placement
controls are `filler_seed` and `filler_min_separation`. `target_conversion` is
in **percent** (`80` means 80%, not `0.8`). Omitting `thickness` produces the bulk
setup; a positive thickness requests the native film geometry.

Crosslink formation and conversion occur during the LAMMPS calculation. The
generator prepares strands, crosslinkers, geometry, and reaction instructions;
it does not claim that the requested conversion has been achieved.

## Coating

```ini
system = coating
model = v22
output_dir = ../runs/coating-PDMS

n1 = 16
m1 = 12
n2 = 16
functionality = 8
m4 = 0
oil = pdms
oil_length = 4
oil_wt = 10
oil_seed = 20260727
density = 0.1
seed = 5489
```

V22 uses the upstream straight-strand formulation; V35 uses its folded-strand
formulation. Selection is explicit because their defaults and reaction settings
differ. `n1` and `m1` are network strand bead length and molecule count; `n2` is
crosslinker bead length. `n4` and `m4` describe five-bead moderators; `m4 = 0`
omits them. Crosslinker count is derived by the backend.

Omit all oil controls for the native no-oil case. Oil-containing cases require
`oil = pdms`, `pmps`, or `copolymer`, plus `oil_length` and `oil_wt`. Copolymer
settings include `mps_percent` or `mps_wt`, `sequence`, and `mps_distribution`.
Weight controls are percentages, and the native integer rounding of molecule
counts is retained. Check the native info record for the achieved composition.

Other settings include `crosslink_distribution`, `crosslink_seed`, `mass`,
`target_density`, `bond_length`, `spacing`, and `oil_min_separation`. A positive
`thickness` creates the film geometry and dependent surface-stage inputs.
Surface controls are `surface_padding`, `surface_relax_steps`,
`surface_production_steps`, and `surface_sample_every`. Omitting thickness
prepares a bulk precursor even when `system = coating`.

## Native settings and provenance

All detailed parameter names, defaults, and validation rules remain those of the
pinned upstream versions. This first integration standardizes selection and
packaging; it does not rename scientific parameters across the three models.
Inspect the complete native help with:

```sh
siliconelab_generator --help oil
siliconelab_generator --help elastomer
siliconelab_generator --help coating-v22
siliconelab_generator --help coating-v35
```

Native help shows `--output` and `--config`; the common interface uses
`output_dir` and a single `--config FILE`. Scientific option names become config
keys, e.g. `--target-density` becomes `target_density`.

For `output_dir = ../runs/my-case`, every system produces:

| File | Purpose |
| --- | --- |
| `data.my-case` | Initial LAMMPS structure and topology. |
| `in.my-case` | Native bulk/network simulation instructions. |
| `submit.my-case.sh` | Inherited cluster submission template. |
| `my-case.info` | Native scientific metadata, resolved defaults, composition, and simulation settings. |
| `request.conf` | Exact copy of the user's common configuration. |
| `generator.conf` | Scientific settings passed to the backend, with an absolute output path. |
| `siliconelab.json` | Schema version, generation status, system/model, upstream repository/commit, and package file inventory. |

Additional film/surface inputs and scripts depend on the system and geometry.
Paths in the manifest's file inventory are relative to the package. The native
info file records actual composition and inherited defaults; `generator.conf`
records supplied settings rather than expanding every default.

Source files are imported without edits and checked against the SHA-256 values
in `upstream.json`. CMake compiles each backend separately, changing only its
entry-point name, and links them into one executable. No network access or source
checkout is needed to build. Original MIT licenses are preserved in `vendor/`.

Generation has been checked against independently compiled native programs for
the included small examples and several advanced parameter combinations. These
checks establish equivalence for those cases, not physical convergence or full
coverage of every scientific configuration. Cluster credentials, job submission,
MD execution, and analysis are subsequent workflow steps.
