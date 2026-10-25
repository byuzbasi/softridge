# Contributing

Use GitHub issues for questions or reproducible bug reports. Include package
versions, a small runnable example and `sessionInfo()`; remove private data
and machine-specific paths. Documentation corrections are welcome.

Changes to numerical formulas, preprocessing, tuning grids or default behavior
must be reviewed explicitly. Preserve compatibility and retain the frozen
reproduction engines and seeds. Use bounded deterministic tests; do not add
production simulations to installation, checks or documentation builds.

Build and check each package locally before proposing a change. Building
the shared site is documented in `tools/build-site.R` and requires pkgdown.
