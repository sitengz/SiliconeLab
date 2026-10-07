# Architecture under evaluation

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

Before implementation, evaluate existing workflow infrastructure for persistent
state, scheduler access, and provenance. Record the reasons for reuse or custom
implementation after examining the three real cases.
