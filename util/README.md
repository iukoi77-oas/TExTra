# Utility Directory Boundary

This directory contains two different kinds of files:

- TExTra project utilities under `util/common/`, `util/prep/`, `util/qual/`, and `util/quant/`
- Vendored external tools under `util/external/`

TExTra project code imports only the project utility packages. The vendored
tool directories are treated as external executables and data bundles:

- TACO is invoked through the `taco_run` executable.
- PLEK2 is invoked through `PLEK2.py` in its own working directory.

Do not import Python modules from `util/external/taco*/` or
`util/external/PLEK*/` into TExTra pipeline code. If project logic is needed,
keep it in the project utility packages or under `TExTra/src/`.

HITindex-related helpers are part of the qual-stage implementation and live
under `util/qual/`.
