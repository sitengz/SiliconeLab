# Contributing to SiliconeLab

The project is starting with three real pilot workflows: silicone oil, PDMS
networks, and silicone coatings. See `workflows/pilots.json` for candidate cases
and `docs/pilot-acceptance.md` for validation expectations.

Use issues to report reproducible problems, discuss scientific requirements, or
propose interfaces. Include the relevant source version, configuration, command,
and a small input/output example when possible. Keep credentials and personal
cluster settings outside submitted examples.

For implementation changes, describe the scientific behavior and its validation
in a pull request. Preserve links to the source tools and make any changed
scientific assumptions explicit. The shared engine and per-system interfaces
will continue to be established through the pilots. The first imported component
is the generator; detailed configuration rules are in `docs/generator.md`.

Generator changes must pass the CMake/CTest checks in the README. Keep vendored
scientific files unchanged unless an intentional scientific change is documented.
When updating an upstream import, record its commit and imported-file checksums
in `upstream.json`, preserve its license, and rerun native parity checks.

## Configurable scientific settings

When changing generators, analyzers, shared headers (`.hpp`), or scientific
workflow logic, expose user-adjustable scientific and analysis choices as named
configuration parameters. Do not tune a shared source variable for one study or
require users to edit and rebuild source to select supported settings. Generator
settings belong in `model.conf`; analyzer settings belong in the corresponding
saved analysis configuration. Existing command-line interfaces may remain, but
new tunable choices must have a reproducible configuration-file interface.

Document each parameter's meaning, units, default, allowed values, applicability
and interactions with other settings. Validate invalid and conflicting inputs,
preserve existing defaults unless a scientific change is intentional, and record
the resolved values (including defaults and derived settings) in run or analysis
metadata. Include a usable configuration example and meaningful tests for
default behavior, an override, and invalid input when adding a parameter.

Keep formulation-specific values in case configurations under `runs/`; keep the
shared implementation reusable. Mathematical constants, algorithm internals and
constraints intrinsic to a supported physical model need not become arbitrary
user controls; document intentional scientific restrictions explicitly. This
policy guides future changes and does not imply every legacy fixed setting
already has a configuration key. Preserve submitted inputs and pinned execution
code when improving configuration interfaces for future runs.

## Repository and run boundaries

Commit reusable generators, analyzers, workflow logic, documentation and tests.
The four material families in `examples/` are shared reference cases for users
to reproduce and test. Keep future experiment configurations, submission state,
Slurm output, trajectories and analysis results under `runs/<batch>/`, which is
ignored by Git. Record the software revision and approved scientific settings
with every batch, and use separate state files to avoid duplicate submissions
or accidentally resuming another batch.

Promote a run configuration to a public example only deliberately, with a clear
scientific purpose, documented output requirements and validation. Remove
credentials and personal cluster state before publishing an example. See the
README for the recommended run directory layout.
