# Changelog

## softridge 0.3.0

- Added opt-in training least-squares rescaling to
  [`hr()`](https://byuzbasi.github.io/softridge/reference/hr.md) and
  [`cv_hr()`](https://byuzbasi.github.io/softridge/reference/cv_hr.md).
  Numerical path, scale and CV score calculations use RcppArmadillo.
- Added fixed positive scale arguments and oracle scale information to
  exact risk and SURE routines. These exact formulas apply to fixed
  tuning and a fixed scale, not a scale estimated from the same
  response.
- Preserved the default unscaled estimator and fold-local threshold
  fractions.
- Retained rank-aware minimum-norm least-squares refitting and
  documented explicit absolute grids for independent implementation
  comparison.

### Shared GitHub documentation (2026-10-25)

- Added a common entry point for softridge and scaleCorrection,
  installation instructions, naming map, worked examples and
  reproduction guidance.
- Added software citations and package-specific reference documentation.
- Preserved numerical functions, package APIs and historical
  reproduction engines.
