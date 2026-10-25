# Soft-thresholding operator

Componentwise soft thresholding, `sign(z) * pmax(abs(z) - gamma, 0)`.

## Usage

``` r
soft_threshold(z, gamma)
```

## Arguments

- z:

  numeric vector or matrix.

- gamma:

  nonnegative threshold.

## Value

An object of the same shape as `z`.

## Examples

``` r
soft_threshold(c(-2, -0.5, 0.3, 1.7), 1)
#> [1] -1.0  0.0  0.0  0.7
```
