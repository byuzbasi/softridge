# Ridge regression path computed through the singular value decomposition

Computes \\(X^T X + k I_p)^{-1} X^T y = X^T (X X^T + k I_n)^{-1} y\\ for
every value in `k`. The singular value decomposition is computed once,
so each additional value of `k` costs \\O(p r)\\ operations with \\r =
rank(X)\\.

## Usage

``` r
ridge_path(X, y, k, sv = NULL)
```

## Arguments

- X:

  prepared (centred/scaled) design matrix.

- y:

  prepared response.

- k:

  positive ridge parameters.

- sv:

  optional precomputed output of the internal SVD helper.

## Value

A `p x length(k)` matrix of ridge coefficients.
