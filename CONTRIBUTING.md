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
