# Repulsive nanoparticle series

This page records the original inspected 32-NP V22 prototype. Its generated
sample files remain intact. New generation uses the shared
[composite engine](composites.md), supporting V22, V35 and generic PDMS matrices,
oil-filled formulations and independently configurable NP interactions. The
legacy config syntax below is still accepted with a new output directory.

The initial prototype accepted `np_count`, `np_radius` (angstrom), and
`np_seed` for an **oil-free V22 bulk** configuration. It was a separate
`coating-nanoparticle` series; the existing oil, elastomer and V22/V35 oil-filled
reference configurations and workflow records retain their original behavior.

```bash
cmake --build build --target siliconelab_generator --parallel 4
./build/siliconelab_generator --config examples/nanoparticles/v22-np32.conf
```

The example prepares `runs/np-repulsive-v22/v22-np32`. A new output directory is
required for each generation. The generator writes initial data and a full
native-length preparation/cure/cooling/equilibration input. It does not run MD
or submit jobs. `workflows/nanoparticles.json` records this series separately.

## Configuration and composition

| Key | Meaning | Default / limits |
|---|---|---|
| `np_count` | Number of rigid FCC particles; activates the extension | Required, 1..512 |
| `np_radius` | Radius used to retain bead centers | 30 angstrom; must contain the first FCC shell and fit the box |
| `np_seed` | NP placement/orientation, precursor packing and NP initialization seed | 20261014; positive signed 32-bit integer |

The NP constituent beads have PDMS-equivalent mass 74 g/mol, zero charge and
atom type **6**. This separate interaction type distinguishes them from neutral
PDMS type 1 without changing their mass or reference LJ parameters. Each
particle has a unique molecule ID following all precursor molecules. NP beads
have no covalent topology and are not reactive crosslinking sites.

At 300 K, sigma=6.44584366 angstrom and epsilon=1.01287845 kcal/mol. FCC
nearest neighbors are at `d=2^(1/6)*sigma=7.2352148776855865` angstrom, with
conventional lattice constant `sqrt(2)*d`. The 30 angstrom-radius bead-center
cut produces 429 beads, matching the separately inspected local template.

The V22 example retains 900 strands of 128 beads, 225 crosslinkers of 32 beads,
and six five-bead moderators: 122,430 precursor beads. Its 32 particles add
13,728 beads for **136,158 total beads and 10.0824042656 wt% NP**. The weight
fraction uses NP mass divided by precursor plus NP mass. These masses are the
agreed PDMS-equivalent prototype masses, not a silica-density calibration.

## Placement and preservation

The initial box volume increases in proportion to total mass so that the native
initial density remains 0.1 g/cm3. NPs are placed first, independently oriented,
with conservative periodic separation. Whole precursor molecules are then
translated and rotated around the particles. Native molecular geometry,
atom IDs, reactive types, bond/angle/dihedral topology, component counts and
stoichiometry are preserved. Wrapped positions have consistent image flags.

NP center separation is at least `2*radius+d`. Precursor beads remain at least
`radius+d` from every NP center under minimum-image distances, which ensures
no initial precursor-NP repulsive contacts. The randomized precursor packing
uses 4.5 angstrom minimum **intermolecular** bead separation; internal bonded
distances remain native. Placement has bounded attempts and fails explicitly
instead of deleting precursor beads or accepting overlaps.

`precursor/` preserves the exact native generated data, input and metadata
before augmentation. The composite `.info` and `nanoparticles.json` record
the expanded box, particle counts, centers, seed and distinct component identity.

## Interactions and rigid integration

NP-NP and NP-precursor interactions act between **constituent beads** as WCA:

```text
U(r) = 4 epsilon [(sigma/r)^12 - (sigma/r)^6] + epsilon, r < 2^(1/6) sigma
U(r) = 0, otherwise
```

NP pair coefficients and cutoffs remain at their 300 K values throughout all
stages, including the native 800 K preparation. For NP pairs with DMS types
1..3 and type 6, the cutoff is 7.2352148776855865 angstrom. The unused MPS types
4/5 have explicit corresponding mixed 300 K WCA coefficients, so all six types
are defined. Internal nonbonded pairs within each NP molecule are excluded.

Precursor pairs retain the exact native 800 K `lj/cut` coefficients/cutoffs and
300 K `lj/gromacs 12 15` coefficients. Native bonded coefficients, crosslinking
probability/seed, timestep, thermostat/barostat targets, compression scale,
native stage lengths and sampling intervals remain unchanged. Three additional
one-step updates (5 fs each) restore rigid coordinates after box changes stop:
one after compression and one before each NPT stage output. The NP pair style is
an additional hybrid substyle; shifting its energy does not shift native pairs.

Each NP translates and rotates through `fix rigid/nvt/small molecule`; ordinary
NVT/NPT integration applies only to the precursor group. A rigid fix is created
before its paired precursor fix, as required by LAMMPS for box changes. Only the
precursor NPT fix controls pressure, using `dilate precursor` for flexible atoms;
LAMMPS separately remaps rigid-body centers. Compression likewise uses the
precursor group. Affine rescaling of NP constituent beads is excluded so that
stage transitions cannot redefine shrunken particles as new rigid bodies.
The rigid state uses `reinit no` within each stage. After box changes stop, one
fixed-volume step refreshes constituent coordinates and image flags before
writing output or redefining the rigid fix. This also handles particles spanning
periodic boundaries; internal geometry is retained across stage transitions.
Its thermostat uses precursor temperature; its pressure compute includes the
entire composite, including NP kinetic and rigid-body virial contributions.
Native affine compression keeps the particles rigid. The communication cutoff
is enlarged to cover each rigid body's owner-to-constituent distance. The added
NPs receive 800 K Gaussian velocities at the recorded seed before first rigid
initialization; native precursor initialization is retained.

This is a **bead-built rigid colloid** model, not the integrated analytical
`pair_style colloid` potential or a single enlarged LJ bead. The cutoff requested
here is the constituent-bead LJ minimum. See the LAMMPS documentation for
[rigid integration](https://docs.lammps.org/fix_rigid.html),
[hybrid pair styles](https://docs.lammps.org/pair_hybrid.html) and
[neighbor exclusions](https://docs.lammps.org/neigh_modify.html).

## Validation and analysis boundary

CTest verifies independent FCC geometry, mass fractions, exact preservation
of native topology and pair/bonded coefficients, covalent distances after
repacking, NP clearance, stage durations, rigid-fix ordering, deterministic
packing, invalid configurations and no-overwrite behavior. Native generator
parity tests remain applicable to configurations without NP settings.

The prepared 32-NP case passed all five CTest suites, an all-stage LAMMPS
preflight, and a 163-step smoke run with LAMMPS 22 Jul 2025 Update 4. All 136,158
beads were retained. Independent audits confirmed exact native topology and
force-field preservation, initial NP clearance and rigid internal distances in
all four stage data outputs. A separate zero-force remapping test placed a
particle across a periodic boundary and verified that the coordinate updates
restore its geometry after compression and NPT. Logs, reproducible audit code
and file checksums are in the case's `validation/` directory. These checks do
not establish convergence or stability of a full production cure.

Two warning categories are recorded in `validation/validation.json`. The NPT
subset-temperature warning is expected: `fix_modify` assigns precursor-only
temperature, then explicitly overrides pressure with the whole-composite
compute. The image-flag warnings arise from newly created type-2 crosslinks
between periodic molecules with different image origins; the audit found no
initial native-bond or NP rigid-body inconsistency. See
[NPT compute overrides](https://docs.lammps.org/fix_nh.html).

The original composite uses the distinct `siliconelab-np-coating-info` metadata
format; new generation uses `siliconelab-composite-info`. The shared structural
analyzer accepts both for full-composite validation and matrix-only basic/topology
analysis. See [composite analysis scope](analyzer.md#composite-basic-and-topology-analysis).
Composite Part 2/Part 3 and NP-obstacle primitive paths remain unsupported;
passing preserved precursor metadata against composite data is invalid.
Generated templates require the LAMMPS
RIGID and native molecular/pair-style packages. Local preflight/smoke outputs
belong under the case's `validation/` directory and are never production results.
