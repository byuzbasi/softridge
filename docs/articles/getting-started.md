# Getting started: two packages, one toolkit

`softridge` is the STR engine; `scaleCorrection` adds the corrected
fitted scale. Both use RcppArmadillo. Install the main package before
its companion.

``` r

remotes::install_github("byuzbasi/softridge", upgrade = "never")
remotes::install_github("byuzbasi/softridge", subdir = "companion/scaleCorrection",
                        upgrade = "never")
```

## Fit STR on a small example

``` r

library(softridge)
set.seed(25)
X <- matrix(rnorm(60 * 100), 60, 100)
beta <- c(2, -1.5, 1, rep(0, 97))
y <- drop(X %*% beta + rnorm(60))
foldid <- rep(1:3, length.out = nrow(X))
grid <- hr_grid(X, y, nk = 3, ng = 5)
fit <- cv_hr(X, y, grid = grid, foldid = foldid)
c(k = fit$k, threshold = fit$gamma, selected = length(fit$selected))
```

    ##            k    threshold     selected 
    ##   0.02851467   0.00000000 100.00000000

``` r

head(fit$coef)
```

    ## [1]  0.01636806  0.89734169 -0.82684181  0.41227578  0.19559541 -0.11110748

The coefficient vector starts with the intercept.
[`hr()`](https://byuzbasi.github.io/softridge/reference/hr.md) returns a
path object; its standard [`coef()`](https://rdrr.io/r/stats/coef.html)
and [`predict()`](https://rdrr.io/r/stats/predict.html) methods select a
path location. The following in-sample predictions illustrate the API,
not prediction accuracy.

``` r

path <- hr(X, y, k = fit$k, gamma = fit$gamma)
head(predict(path, X))
```

    ## [1]  2.5790156 -0.5012138 -3.2782377  2.5351456 -1.8415457  0.8417265

## Manuscript names and code entry points

| Manuscript name | Code | Scale |
|:---|:---|:---|
| STR | [`hr()`](https://byuzbasi.github.io/softridge/reference/hr.md) / [`cv_hr()`](https://byuzbasi.github.io/softridge/reference/cv_hr.md) | None |
| STR-OLS | `cv_hr(..., refit = TRUE)` | Selected-variable LS refit |
| SSTR-LS | `cv_hr(..., scale = TRUE)` | Training LS scale |
| SSTR | `scaleCorrection::corrected_pair()$SSTR_c0` | Corrected scale |
| Ridge-c0 | `scaleCorrection::corrected_pair()$Ridge_c0` | Corrected zero-threshold control |

## Small-example provenance

This 60-by-100 Gaussian example uses seed 25, three supplied folds and a
3-by-5 grid. It demonstrates software interfaces. It is not one of the
reported performance comparisons and does not establish a method
ranking.

## Threshold geometry

The illustration below evaluates the deterministic soft-threshold
operator; it contains no estimated performance results. Source:
`tools/build-gallery.R`.

![Piecewise linear soft threshold operator at threshold
one](../reference/figures/threshold-geometry.png)
