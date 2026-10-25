# Default tuning grids for the STR estimator

The ridge grid is geometric between `kmin_ratio * d1^2` and `d1^2`,
where `d1` is the largest singular value of the prepared design. For
each ridge value the threshold grid is `c * max|ridge|`, `c` in
`[0, 0.99]`.

## Usage

``` r
hr_grid(
  X,
  y,
  nk = 12,
  ng = 40,
  kmin_ratio = 1e-04,
  intercept = TRUE,
  standardize = TRUE
)
```

## Arguments

- X:

  design matrix.

- y:

  response.

- nk:

  number of ridge values.

- ng:

  number of threshold values.

- kmin_ratio:

  smallest ridge value relative to the squared leading singular value.

- intercept, standardize:

  passed to the design preparation.

## Value

A list with the vector `k` and the matrix `gamma` (`nk x ng`), one row
of thresholds per ridge value, and `fraction`, recomputed on the
training response inside each CV fold.
