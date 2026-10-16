# SiliconeLab development instructions

Follow the configurable scientific settings policy in [CONTRIBUTING.md](CONTRIBUTING.md).
For generator, analyzer, shared-header and scientific workflow changes, expose
user-adjustable choices through saved configuration parameters instead of
editing shared constants for a particular study. Document and validate the
parameters, preserve established defaults, record resolved values in output
metadata, and verify default and override behavior. Keep experiment-specific
settings under `runs/` and preserve submitted inputs and pinned executors.
