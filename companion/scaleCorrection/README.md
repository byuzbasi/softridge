# scaleCorrection

**Optimism-corrected scaling for soft-thresholded ridge.**

This RcppArmadillo companion depends on the matching `softridge` engine.
It jointly tunes corrected SSTR and the zero-threshold Ridge-c0 control
using the same observations, outer folds and ridge-penalty candidates.

Install `softridge` first, then install this subdirectory:

```r
remotes::install_github("byuzbasi/softridge", subdir = "companion/scaleCorrection",
                        upgrade = "never")
```

See the [shared getting-started guide](https://byuzbasi.github.io/softridge/articles/getting-started.html),
[worked scaling example](https://byuzbasi.github.io/softridge/articles/scaling.html),
and [function reference](https://byuzbasi.github.io/softridge/scaleCorrection/reference/index.html).

## Corrected scale and matched control

For a centered training response and the fitted vector from thresholded ridge,
the scale is the nonnegative part of the response-fit inner product minus
estimated variance times the selected ridge diagonal sum, divided by the
squared norm of the fitted vector. A zero fitted vector receives scale zero.
The returned `SSTR_c0` and `Ridge_c0` records contain intercept-first
coefficient vectors, selected penalties, scales and CV scores.

The variance estimate uses training-only ten-fold CV-Lasso and records its
fold digest. This plug-in correction and response-adaptive tuning do not
inherit a fixed-parameter unbiased-risk guarantee.

| Return component | Meaning |
|:--|:--|
| `SSTR_c0` | Jointly tuned ridge penalty and threshold with corrected scaling |
| `Ridge_c0` | Zero-threshold control tuned over the same ridge penalties |
| `cv_error` | Mean validation squared-error surface |
| `variance` | Training-subset variance diagnostics and fold digests |
| `foldid`, `grid` | User-supplied outer folds and candidate definitions |

The methods and fixed data-processing rules are those of the frozen research
implementation, version 0.1.0. This package is distributed on GitHub, not CRAN.
