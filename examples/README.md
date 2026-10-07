# Generator examples

These small configurations exercise preparation and are not production-size
scientific pilots:

| Config | System |
| --- | --- |
| `oil.conf` | PDMS oil, N16, 12 chains. |
| `elastomer.conf` | Linear PDMS strands with tetrafunctional crosslinkers. |
| `coating-v22.conf` | V22 straight-strand network with PDMS oil. |
| `coating-v35.conf` | V35 folded-strand network with copolymer oil. |

From the repository root, run
`./build/siliconelab_generator --config examples/elastomer.conf`, for example.
Outputs go to `runs/` beside the repository's `examples/` directory. Each output
directory must be new. Edit `output_dir` for another run.

See `../docs/generator.md` for configuration details. Larger candidate production
pilots remain recorded in `../workflows/pilots.json` and require separate
simulation and analysis validation.

The four [full-size follow-up configurations](full-size/README.md) reproduce
existing production cases after the small Nova jobs pass their completion checks.
