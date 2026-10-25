# Corrected SSTR and the matched Ridge control

Corrected SSTR applies a nonnegative optimism correction to the LS
scale.
[`corrected_pair()`](https://byuzbasi.github.io/softridge/scaleCorrection//reference/corrected_pair.html)
selects a thresholded fit and the zero-threshold Ridge control on the
same outer folds and ridge penalties. This control isolates the
empirical contribution of thresholding within the scaled family.

``` r

library(softridge)
library(scaleCorrection)
set.seed(25)
X <- matrix(rnorm(60 * 100), 60, 100)
y <- drop(X[, 1:3] %*% c(2, -1.5, 1) + rnorm(60))
grid <- hr_grid(X, y, nk = 3, ng = 5)
foldid <- rep(1:3, length.out = nrow(X))
pair <- corrected_pair(X, y,
  grid = list(k = grid$k, fraction = grid$fraction),
  foldid = foldid, unit_id = "showcase-25")
pair$SSTR_c0[c("k", "gamma", "scale", "size")]
pair$Ridge_c0[c("k", "gamma", "scale", "size")]
prediction <- pair$SSTR_c0$coef[1] + X %*% pair$SSTR_c0$coef[-1]
```

The same example is executable as `tools/smoke-companion.R` from the
repository root after installing both packages. It is intentionally
bounded. Its variance estimator uses ten-fold CV-Lasso in each supplied
training subset; this is additional to the three outer folds in this
demonstration.

## Returned quantities

`SSTR_c0` and `Ridge_c0` contain intercept-first `coef`, ridge penalty
`k`, threshold `gamma`, scale, selected size and validation score.
`cv_error` stores the shared mean-squared validation surface. `variance`
records training-only variance estimates and fold hashes. Supplied
`foldid` and `grid` are returned for traceability.

The correction uses the selected ridge diagonal sum as a plug-in
optimism term. It is a practical extension of established scaling ideas,
including Hagiwara’s orthogonal-model work. An adaptive scale and
selected tuning parameters do not inherit the fixed-tuning unbiased-risk
statement. See [risk and
targets](https://byuzbasi.github.io/softridge/articles/risk-and-targets.md)
for the separate quantities.
