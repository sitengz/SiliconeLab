# SiliconeLab

A scientific workflow project for silicone oils, PDMS networks, and silicone coatings.

SiliconeLab provides a C++17 generator selected by `system = oil`, `elastomer`,
or `coating` in a saved configuration. It reuses the original scientific
generators and writes a common run package with configuration and source records.

The structural analyzer accepts an equilibrated data snapshot and its matching
info file for oil, elastomer and V22/V35 coating systems. Generalized part 2
and part 3 analyzers cover spatial/dynamic and surface/mechanical properties,
with explicit applicability and sampling checks. The four reference cases
exercise the shared tools; future simulation batches belong under `runs/`.

## Generate a system

Build with a C++17 compiler and CMake 3.19 or newer:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build --config Release --parallel
./build/siliconelab_generator --config examples/oil.conf
```

On a Windows multi-configuration build, the executable is
`build/Release/siliconelab_generator.exe`. Generation requires no Python or AI
service. The CMake build is intended for Linux, macOS, and Windows; WSL is the
initial locally verified environment.

```ini
system = oil
output_dir = ../runs/my-oil

length = 16
chains = 12
mps_percent = 0
seed = 20260727
```

`output_dir` resolves relative to the configuration file and must be new. Detailed
settings use the selected generator's existing keys and defaults. Coating also
requires `model = v22` or `model = v35`. See [the configuration guide](docs/generator.md)
and [examples](examples/README.md) for all three systems.

The output includes a LAMMPS data file, simulation inputs, inherited Slurm
templates, native scientific metadata, the original and backend configs, and
`siliconelab.json`. Generation does not submit jobs or run molecular dynamics.
Slurm templates still contain the upstream cluster assumptions and need review
before use on another cluster.

## Analyze a system

```sh
./build/siliconelab_analyzer data.CASE.npt_eq CASE.info --z1-command /path/to/Z1+
```

Twelve structural/network analyses are supported. C++17 performs geometry and
graph calculations; Python 3.9+ coordinates validation, reports and external Z1+.
The analyzer writes a new output directory and preserves simulation files.
Network quantities are marked not applicable for uncrosslinked oil. See
[the analysis guide](docs/analyzer.md) for definitions, component selection,
Z1 contour mapping, outputs and boundary limitations.

## Initial scope

The first milestone is one validated workflow from each source repository:

| System | Source | Candidate pilot | Requested analysis |
| --- | --- | --- | --- |
| Oil | sitengz/Silicone_Oil | N16_PDI1 | Surface tension and sampling diagnostics |
| Network | sitengz/PDMS_Elastomer | Existing N40 linear-network case | Conversion and surface profiles |
| Coating | sitengz/Silicone_Coating | Existing V22 oil-containing case | Select a property supported by its current analyzer |

The four small examples and four full-size bulk configurations have completed
on Nova through private local orchestration profiles. A portable persistent
submission and monitoring controller remains a subsequent milestone. Matched films and dependent surface, dynamics and tensile stages are defined in
[the film workflow](workflows/films.json), with persistent accepted-job records
and output gates. See [parts 2 and 3](docs/analysis-parts-2-3.md) for analysis inputs
and the two generalized analyzer commands.

## Planned workflow

Scientific instruction -> explicit job specification -> validation -> generation
-> cluster submission -> monitoring -> completion checks -> output selection
-> retrieval and analysis -> figures, tables, diagnostics, and provenance.

Molecular dynamics remains a required calculation. Completion of a job does not
establish equilibration or convergence.

## Layout

- `workflows/`: published reference workflow specifications and pilot records.
- `profiles/`: example cluster settings; real credentials stay outside the repository.
- `docs/`: architecture, operating permissions, and pilot acceptance criteria.
- `src/`: generalized C++ generator/analysis backends and Python analysis controllers.
- `scripts/`: preparation, submission, verification and analysis orchestration.
- `vendor/`: unchanged scientific sources with their original MIT licenses.
- `upstream.json`: pinned source versions and imported-file checksums.
- `engine/`: reserved for persistent workflow execution and job records.
- `modules/`: scientific integration notes.
- `examples/`: four reference cases: oil, elastomer, V22 and V35, with small,
  full-size bulk and matched-film configurations.
- `runs/`: ignored, batch-specific configurations, inputs, job records, trajectories
  and analysis reports for future simulations.
- `tests/`: scientific contract and parity checks against native programs.

The existing repositories remain the scientific implementation sources. Imported
generator files retain their contents, attribution, licenses, and links to their
original commit history.

## Shared software and future runs

GitHub contains the reusable generator, all three analyzer entry points,
workflow scripts, documentation, tests and four reference cases. The four cases
are reproducible examples for other users to test the tools, rather than a limit
on the configurations the tools support. Large simulation outputs, credentials
and machine-specific state are excluded from Git.

Keep reference configurations in `examples/` stable. For a future scientific
batch, copy the relevant example into a new directory under `runs/` and record
its approved settings there. Use a unique batch and case name:

```text
runs/<batch>/
  workflow.json                 # approved cases, stages, resources and outputs
  workflow-state.json           # accepted job IDs and monitoring state
  <case>/
    configs/model.conf          # complete approved system configuration
    generated/                  # generator package and provenance
    surface/                    # stage inputs, Slurm output and simulation output
    dynamics/
    tensile-x/                  # only when applicable
    tensile-y/
    analysis/part1/
    analysis/part2/
    analysis/part3/
```

This is the directory convention for new batches, not a new command-line API.
The existing reference film controller currently uses `workflows/films.json`
and `runs/film-workflow.json`; active reference jobs retain those paths. A future
batch must use its own plan and state paths before submission so it cannot
reuse the reference batch's job records.

The generator resolves `output_dir` relative to the configuration file and
requires a new output directory. For the layout above, use
`output_dir = ../generated` inside `configs/model.conf`. Preserve the submitted
configuration, code revision, seeds, checksums, dependency records and analysis
requirements with each run. Submit from each stage's directory and keep its
Slurm stdout/stderr there. On Nova, dated activity summaries remain in the
central log directory.

Only promote a future case into `examples/` deliberately, after documenting
its purpose and validation. Changes to shared scientific code, workflow logic,
analysis or documentation still belong in GitHub. See
[the example guide](examples/README.md) for copying a reference configuration.

## Development status

Local generation and structural, spatial/dynamic, surface and tensile analysis
controllers are implemented. Film properties retain explicit applicability and
sampling/geometry quality checks. Nova pilots
have completed four small and four full-size bulk simulations using private
local orchestration profiles. Retrieval, a portable persistent controller,
recovery and the AI interface remain future work.

To run the generator and analyzer checks (Python 3.9+ is required for testing):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```
