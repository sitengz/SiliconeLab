# Scientific modules

The generator now selects `oil`, `elastomer`, or `coating` using the unchanged
source implementations under `../vendor/`. The elastomer backend comes from
PDMS_Elastomer, recorded as `network` in `../upstream.json`. Coating supports both
V22 and V35. Common selection and packaging live in `../src/`.

The generalized static analyzer maps oil, elastomer and coating metadata onto
the pinned C++ structural calculations, with component validation and external
Z1+ integration. See `../docs/analyzer.md`. Persistent portable simulation
orchestration remains future work; machine-specific Nova pilots use local profiles.
