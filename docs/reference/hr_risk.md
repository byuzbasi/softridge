# Exact finite-sample moments and risk of the STR estimator

Under \\Y = X\beta + \epsilon\\, \\\epsilon \sim N_n(0, \sigma^2 I_n)\\
and a fixed design, returns the exact mean, second moment and mean
squared error of each coordinate of the STR estimator as an estimator of
the projection \\\theta = P_X \beta\\, and the total risk
\\E\\\hat\theta - \theta\\\_2^2\\. No centring or scaling is applied:
`X` is used as supplied.

## Usage

``` r
hr_risk(X, beta, sigma, k, gamma, nonneg = FALSE, scale = 1)
```

## Arguments

- X:

  design matrix.

- beta:

  true coefficient vector.

- sigma:

  error standard deviation.

- k:

  ridge parameter.

- gamma:

  threshold.

- nonneg:

  logical, nonnegative variant.

- scale:

  For `hr`, logical training least-squares rescaling (SSTR-LS). For
  `hr_risk`, a fixed positive numeric scale; not a random fitted scale.

## Value

A list with components `mean`, `second`, `mse`, `risk`, `theta`, `m`
(mean of the ridge estimator) and `s` (its coordinatewise standard
deviations).
