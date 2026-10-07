# Architecture

The first implemented component is `siliconelab_generator`, a C++17 executable
with a common `system`, coating `model`, and `output_dir` config interface. Each
scientific backend is compiled separately from unchanged, pinned upstream source
and linked in-process. The dispatcher removes common keys and passes detailed
settings through a saved native config. See `generator.md` for the contract.

All systems write to the same package convention regardless of their original
output directory rules. The package preserves requested settings, native resolved
metadata, backend source version, and generated file inventory. Existing packages
are rejected; generation failures remove only the package created by that call.

This component makes no scheduler calls. Generated submission scripts retain the
upstream cluster settings until cluster-profile integration is implemented.

The common engine will store run state and coordinate preparation, submission,
monitoring, completion checks, output selection, transfer, and analysis.

Oil, network, and coating modules will define their scientific settings,
generator commands, simulation stages, expected outputs, and analyzer commands.
The pilots determine which interfaces can be shared and which remain specific.

The optional AI interface produces an explicit, reviewable job specification.
Execution must also be possible from saved specifications without an AI service.

Run records should link the original instruction, interpreted specification,
configurations, source commits, seeds, input files, job IDs, stage outcomes,
selected outputs, analyzer versions, and final results.

Cluster profiles are machine-specific. Credentials are provided by the local
authentication mechanism and are not stored in specifications or source control.

Before implementing cluster orchestration, evaluate existing infrastructure for persistent
state, scheduler access, and provenance. Record the reasons for reuse or custom
implementation after examining the three real cases.
