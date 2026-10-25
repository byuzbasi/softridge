# Validation scope

The public source preserves the frozen research R functions, C++
kernels, default parameters and input preprocessing. Documentation and
repository metadata are additions. Source fingerprints are recorded
locally before publication; the original versioned tarballs remain in
`reproduction/packages/`.

The supplied manuscript validation compared the RcppArmadillo
implementations with the reference formulas on two saved matrix
fixtures: 107 checks, maximum absolute difference 8.65e-12. This is
checked-formula agreement, not a claim that the reference and
independent studies use identical default CV protocols.

Publication checks include local package build/check, the original
softridge test suite, a bounded corrected-scale example, source/checksum
validation, portable reproduction help/checks, documentation rendering
and link/privacy inspection. Exact check outcomes are listed below after
validation.

No full numerical study is launched to build this repository or its
website. No CRAN publication, Windows check, cloud CI result or speedup
is implied. The scientific claims concern fixed-parameter formulas;
practical scaling and cross-validation are evaluated separately in the
companion manuscript.
