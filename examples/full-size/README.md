# Full-size workflow cases

These four configurations reproduce production cases from the source repositories.
They follow the four small generator examples after those jobs finish successfully.

| Configuration | Production system |
| --- | --- |
| `oil.conf` | PDMS oil, N16, 6,250 chains: 100,000 beads. |
| `elastomer.conf` | Linear N40 network, 3,000 strands, N32 tetrafunctional crosslinkers, 95% target conversion: 168,000 beads. |
| `coating-v22.conf` | V22 N128 network, 900 strands, octafunctional crosslinkers, 10 wt% PDMS oil with N30. |
| `coating-v35.conf` | V35 N384 folded network, 306 strands, tetrafunctional crosslinkers, 10 wt% random copolymer oil with N15 and 50 mol% MPS. |

The coating cases retain the source formulation's six five-bead moderators.
The requested oil weight fraction is realized with whole oil molecules; exact
counts and realized composition are recorded by the generator. The elastomer
retains the existing conversion-controlled protocol rather than the small
example's fixed-duration crosslinking protocol.

Source configuration paths, checksums, pinned generator versions, required
outputs, and Nova resource settings are recorded in
[`workflows/full-size.json`](../../workflows/full-size.json).

The launch sequence is:

1. Verify that all four small jobs completed with exit `0:0`, completion markers,
   final equilibrated data, and the expected network/coating MSD trajectories.
2. Commit these configurations and the workflow record to the existing development
   branch, push to GitHub, and pull that exact revision into the Nova checkout.
3. Generate the full-size packages in a short preparation job on a compute node.
   Check its exit status and all four package manifests before simulation submission.
4. Submit four full simulations from their respective generated run folders.
   Each requests one node, 96 MPI ranks, 200 GB, and 48 hours, matching the prior
   Nova templates. Keep Slurm stdout/stderr in that case folder.
5. Monitor hourly using Slurm accounting and output checks. Keep independent jobs
   within the approved limit of eight and prevent duplicate submissions.

Generation uses the standard interface, for example:

```sh
./build/siliconelab_generator --config examples/full-size/oil.conf
```

Outputs go to separate new directories under `runs/`; completed small runs are
preserved. Full-size generation, compilation, analysis, and molecular dynamics
run on compute nodes. Preserve the full generated scientific inputs and stage
durations. This batch contains the four primary bulk workflows; oil-film and
coating surface stages are separate future runs.

A job finishing successfully demonstrates execution and output continuity.
Equilibration, conversion attainment, and property convergence require scientific
checks after the outputs are available.
