# SiliconeLab

A scientific workflow project for silicone oils, PDMS networks, and silicone coatings.

SiliconeLab provides a C++17 generator selected by `system = oil`, `elastomer`,
or `coating` in a saved configuration. It reuses the original scientific
generators and writes a common run package with configuration and source records.

Cluster simulations, output selection, property analysis, and an optional AI
interface are subsequent milestones.

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

## Initial scope

The first milestone is one validated workflow from each source repository:

| System | Source | Candidate pilot | Requested analysis |
| --- | --- | --- | --- |
| Oil | sitengz/Silicone_Oil | N16_PDI1 | Surface tension and sampling diagnostics |
| Network | sitengz/PDMS_Elastomer | Existing N40 linear-network case | Conversion and surface profiles |
| Coating | sitengz/Silicone_Coating | Existing V22 oil-containing case | Select a property supported by its current analyzer |

Pilot configurations and analysis settings must be checked against the source
tools before execution. No cluster connection or simulation submission is
implemented yet. The small generator examples are preparation checks, separate
from these production pilots.

## Planned workflow

Scientific instruction -> explicit job specification -> validation -> generation
-> cluster submission -> monitoring -> completion checks -> output selection
-> retrieval and analysis -> figures, tables, diagnostics, and provenance.

Molecular dynamics remains a required calculation. Completion of a job does not
establish equilibration or convergence.

## Layout

- `workflows/`: system-specific workflow specifications and pilot records.
- `profiles/`: example cluster settings; real credentials stay outside the repository.
- `docs/`: architecture, operating permissions, and pilot acceptance criteria.
- `src/`: the common C++ generator entry point and backend adapters.
- `vendor/`: unchanged generator sources with their original MIT licenses.
- `upstream.json`: pinned source versions and imported-file checksums.
- `engine/`: reserved for persistent workflow execution and job records.
- `modules/`: scientific integration notes.
- `examples/`: small configurations for oil, elastomer, V22, and V35.
- `tests/`: parity checks against independently compiled native generators.

The existing repositories remain the scientific implementation sources. Imported
generator files retain their contents, attribution, licenses, and links to their
original commit history.

## Development status

Local generation is implemented. Cluster authentication, submission, monitoring,
retrieval, analyzers, recovery, and the AI interface remain to be established
through the three pilots.

To run the generator checks (Python 3 is required for testing):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```
