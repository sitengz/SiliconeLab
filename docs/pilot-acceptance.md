# Pilot acceptance

For each pilot, record the scientific question, exact configuration, source code
version, build and execution commands, seeds, cluster settings, required outputs,
selection rules, and analysis settings.

1. Prepare a manual reference workflow and verify the expected inputs and outputs.
2. Reproduce generation with the same configuration and seed where supported.
3. Run the actual required simulation stages through the scheduler.
4. Record job IDs and stage dependencies; prevent duplicate submissions.
5. Check scheduler outcomes, simulation logs, required outputs, and stage continuity.
6. Verify selected outputs and transfers before analysis.
7. Compare analysis of identical saved data against the manual reference using
   numerical tolerances appropriate to the calculation.
8. Assess independent trajectories with an appropriate statistical comparison;
   do not require trajectories to be numerically identical.
9. Report scientific quality concerns separately from execution success.
10. Demonstrate that controller restart preserves job state and that a failed
    stage is reported without silently changing scientific settings.

Short runs may validate execution plumbing. Scientific estimates need the actual
sampling protocol and quality checks; a short smoke run is not production validation.

Track human interventions and delays between stages. Queue time and molecular
dynamics runtime should be reported separately from automation overhead.
