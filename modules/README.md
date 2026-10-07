# Scientific modules

The generator now selects `oil`, `elastomer`, or `coating` using the unchanged
source implementations under `../vendor/`. The elastomer backend comes from
PDMS_Elastomer, recorded as `network` in `../upstream.json`. Coating supports both
V22 and V35. Common selection and packaging live in `../src/`.

Simulation orchestration and analyzer integration remain future work.
