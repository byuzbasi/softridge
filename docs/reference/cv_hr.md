# Cross-validated STR estimator

Selects the ridge parameter and the threshold by K-fold cross-validation
of the prediction error. With `refit = TRUE`, the variables retained by
thresholding are refitted by least squares and the pair (ridge
parameter, model size) is selected instead; model size `m` corresponds
to the data-dependent threshold equal to the `(m+1)`-th largest absolute
ridge coefficient when magnitudes are distinct. Ties use column order.
Refit size is capped by the smallest training-fold rank, and refits use
the Moore-Penrose solution on the prepared design (relative tolerance
1e-10).

## Usage

``` r
cv_hr(
  X,
  y,
  nfolds = 10,
  grid = NULL,
  refit = FALSE,
  max_size = NULL,
  foldid = NULL,
  intercept = TRUE,
  standardize = TRUE,
  scale = FALSE
)
```

## Arguments

- X:

  design matrix.

- y:

  response.

- nfolds:

  number of folds.

- grid:

  output of `hr_grid`; computed on the full data if missing.

- refit:

  logical, least-squares refit on the selected variables.

- max_size:

  largest model size considered when `refit = TRUE`.

- foldid:

  optional fold assignment.

- intercept, standardize:

  passed to the design preparation.

- scale:

  Logical, estimate a least-squares scale inside each training fold and
  on the final training sample (SSTR-LS). Incompatible with refit.

## Value

A list with the selected tuning values, the cross-validation error
surface, the final coefficient vector `coef` (intercept first) and the
selected index set `selected`.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(60 * 150), 60, 150)
y <- X[, 1:3] %*% c(3, -2, 1.5) + rnorm(60)
cv <- cv_hr(X, y, nfolds = 5)
cv$selected
#>  [1]   1   2   3   7  11  14  16  19  21  22  24  30  32  33  36  38  39  40  42
#> [20]  43  45  47  49  54  56  57  59  61  64  65  66  67  70  71  72  76  77  81
#> [39]  82  83  95 101 105 110 111 112 114 115 116 117 118 119 121 123 126 128 130
#> [58] 131 133 134 135 137 139 140 142 145 149
```
