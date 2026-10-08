# Reference cases

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

The [matched-film configurations](films/) use the same four material families.
Bulk and film examples document reference settings and analysis requirements;
cluster resources and historical prerequisite job IDs in reference workflow
plans must be reviewed for another environment.

## Future simulation configurations

Copy a reference configuration into `runs/<batch>/<case>/configs/model.conf`
and make scientific changes in that copy. Keep the shared examples intact.
For instance, from the repository root:

```sh
mkdir -p runs/my-batch/my-oil/configs
cp examples/oil.conf runs/my-batch/my-oil/configs/model.conf
```

Set `output_dir = ../generated` in the copied configuration, then generate:

```sh
./build/siliconelab_generator --config runs/my-batch/my-oil/configs/model.conf
```

The generated package goes to `runs/my-batch/my-oil/generated/`. That directory
must not already exist. Store the batch's workflow plan, accepted job records,
outputs and reports under the same batch directory. `runs/` is excluded from
GitHub; reusable code and deliberately published reference cases remain tracked.
