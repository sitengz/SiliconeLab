# SiliconeLab

A scientific workflow project for silicone oils, PDMS networks, and silicone coatings.

SiliconeLab is being prepared to connect configuration-driven C++ generators,
cluster simulations, output selection, and property analysis through one recorded
workflow. An optional AI interface will help translate scientific questions into
explicit configurations and analysis plans.

## Initial scope

The first milestone is one validated workflow from each source repository:

| System | Source | Candidate pilot | Requested analysis |
| --- | --- | --- | --- |
| Oil | sitengz/Silicone_Oil | N16_PDI1 | Surface tension and sampling diagnostics |
| Network | sitengz/PDMS_Elastomer | Existing N40 linear-network case | Conversion and surface profiles |
| Coating | sitengz/Silicone_Coating | Existing V22 oil-containing case | Select a property supported by its current analyzer |

Pilot configurations and analysis settings must be checked against the source
tools before execution. No cluster connection or simulation submission is
implemented in this scaffold.

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
- `upstream.json`: source repository versions and integration status.
- `engine/`: reserved for persistent workflow execution and job records.
- `modules/`: reserved for oil, network, and coating integrations.
- `examples/`: reserved for reproducible user examples.

The existing repositories remain the scientific implementation sources while
their workflows are evaluated. Consolidation will preserve attribution and
development history, with source licenses checked before importing code.

## Development status

This is a project scaffold, not an operational agent. Execution interfaces,
cluster authentication, resource limits, recovery rules, and quality criteria
will be established through the three pilots.
