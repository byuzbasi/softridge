# Stein unbiased risk estimate and degrees of freedom of the STR fit

For a fixed design without intercept, returns the degrees of freedom
\\\sum\_{j \in \hat S} (A_k)\_{jj}\\, with \\A_k = (X^T X + k I)^{-1}
X^T X\\, and the unbiased estimate of the prediction risk
\\E\\X\hat\theta - X\beta\\^2\\ for fixed tuning and known error
variance. Plug-in variance estimates or tuning by minimizing SURE do not
retain that guarantee. When `intercept = TRUE` the design and response
are centred and one degree of freedom is added for the intercept.

## Usage

``` r
hr_sure(
  X,
  y,
  sigma2,
  k,
  gamma,
  intercept = TRUE,
  standardize = TRUE,
  scale = 1
)
```

## Arguments

- X:

  design matrix.

- y:

  response.

- sigma2:

  error variance (known or estimated).

- k:

  ridge parameters.

- gamma:

  thresholds.

- intercept:

  logical.

- standardize:

  logical.

- scale:

  Fixed positive numeric scale of the slope fit. Intercept df remains
  one. These formulas do not apply to a scale estimated from the same
  response.

## Value

A list with matrices `sure` and `df` of dimension
`length(k) x length(gamma)`.
