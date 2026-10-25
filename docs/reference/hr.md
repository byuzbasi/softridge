# Soft-thresholded ridge (STR) estimator

Fits the soft-thresholded ridge estimator \\\hat\theta_j =
sign(\hat\beta\_{k,j})(\|\hat\beta\_{k,j}\| - \gamma)\_+\\, where
\\\hat\beta_k\\ is the ridge estimator with ridge parameter `k`. The
estimator has an explicit form and is valid for \\p \> n\\.

## Usage

``` r
hr(
  X,
  y,
  k,
  gamma,
  intercept = TRUE,
  standardize = TRUE,
  nonneg = FALSE,
  scale = FALSE
)

# S3 method for class 'hr'
coef(object, ik = 1, ig = 1, ...)

# S3 method for class 'hr'
predict(object, newx, ik = 1, ig = 1, ...)
```

## Arguments

- X:

  `n x p` design matrix.

- y:

  response vector of length `n`.

- k:

  ridge parameter(s), positive.

- gamma:

  threshold(s), nonnegative. Thresholds apply on the scale of the
  standardized design when `standardize = TRUE`.

- intercept:

  logical, fit an unpenalized intercept by centring.

- standardize:

  logical, scale columns of `X` to unit variance.

- nonneg:

  logical, use the nonnegative variant \\(\hat\beta\_{k,j} -
  \gamma)\_+\\.

- scale:

  For `hr`, logical training least-squares rescaling (SSTR-LS). For
  `hr_risk`, a fixed positive numeric scale; not a random fitted scale.

- object:

  a fitted `"hr"` object.

- ik, ig:

  indices of the ridge parameter and threshold to use.

- ...:

  unused.

- newx:

  matrix of new predictor values.

## Value

An object of class `"hr"` containing the coefficient array `beta` of
dimension `p x length(k) x length(gamma)` on the original scale, the
intercepts `a0`, and the tuning grids.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(50 * 200), 50, 200)
y <- X[, 1:3] %*% c(3, -2, 1.5) + rnorm(50)
fit <- hr(X, y, k = 10, gamma = 0.1)
which(coef(fit)[-1] != 0)
#>  [1]   1   2   3  10  11  18  19  22  25  26  27  28  37  39  44  47  49  52  55
#> [20]  60  65  66  68  73  78  81  88  90  92  96  97  98 100 101 102 105 113 115
#> [39] 119 121 122 126 128 130 132 134 136 138 140 144 148 149 151 162 163 164 165
#> [58] 167 172 174 175 176 181 182 185 191 192 193 197 199 200
```
