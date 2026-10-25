# Validation scope

The public source preserves the frozen research R functions, C++ kernels,
default parameters and input preprocessing. Documentation and repository
metadata are additions. Source fingerprints are recorded locally before
publication; the original versioned tarballs remain in `reproduction/packages/`.

The supplied manuscript validation compared the RcppArmadillo implementations
with the reference formulas on two saved matrix fixtures: 107 checks, maximum
absolute difference 8.65e-12. This is checked-formula agreement, not a claim
that the reference and independent studies use identical default CV protocols.

Publication checks include local package build/check, the original softridge
test suite, a bounded corrected-scale example, source/checksum validation,
portable reproduction help/checks, documentation rendering and link/privacy
inspection. Exact check outcomes are listed below after validation.

No full numerical study is launched to build this repository or its website.
No CRAN publication, Windows check, cloud CI result or speedup is implied.
The scientific claims concern fixed-parameter formulas; practical scaling
and cross-validation are evaluated separately in the companion manuscript.

## Local checks for this publication (2026-10-25)

* Both source packages built, installed and passed local `R CMD check
  --no-manual` with status OK (zero errors, warnings and notes).
* The original softridge tests and executable vignette examples passed.
* A single 60-by-100 corrected SSTR/Ridge-c0 example passed its dimension,
  finite-score, nonnegative-scale and zero-threshold-control checks.
* All R expressions and C++ kernels match the frozen research implementation;
  documentation and repository metadata are the only package-source edits.
* No package-name collision was found in current/archived CRAN (2026-10-25)
  or the checked Bioconductor release indexes (2026-10-24). This is an
  availability check, not a CRAN submission or name reservation.

Environment: R 4.6.0 on macOS ARM64; Rcpp 1.1.1.1.1,
RcppArmadillo 15.2.7.1, glmnet 5.0 and digest 0.6.39.
