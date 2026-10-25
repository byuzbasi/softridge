# STR estimator tuned by the Stein unbiased risk estimate

Minimizes the plug-in SURE criterion over the grid. Its selected minimum
is not an unbiased estimate of the risk of the tuned estimator.

## Usage

``` r
sure_hr(X, y, sigma2, grid = NULL, intercept = TRUE, standardize = TRUE)
```

## Arguments

- X:

  design matrix.

- y:

  response.

- sigma2:

  error variance or an estimate of it.

- grid:

  output of `hr_grid`; computed if missing.

- intercept, standardize:

  passed to the design preparation.

## Value

A list with the selected tuning values, the fitted coefficients
(intercept first), the selected set and the degrees of freedom.
