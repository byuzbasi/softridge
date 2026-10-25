# Package index

## Fit STR and choose tuning parameters

- [`hr()`](https://byuzbasi.github.io/softridge/reference/hr.md)
  [`coef(`*`<hr>`*`)`](https://byuzbasi.github.io/softridge/reference/hr.md)
  [`predict(`*`<hr>`*`)`](https://byuzbasi.github.io/softridge/reference/hr.md)
  : Soft-thresholded ridge (STR) estimator
- [`cv_hr()`](https://byuzbasi.github.io/softridge/reference/cv_hr.md) :
  Cross-validated STR estimator
- [`hr_grid()`](https://byuzbasi.github.io/softridge/reference/hr_grid.md)
  : Default tuning grids for the STR estimator

## Fixed-parameter risk and SURE

- [`hr_risk()`](https://byuzbasi.github.io/softridge/reference/hr_risk.md)
  : Exact finite-sample moments and risk of the STR estimator
- [`hr_sure()`](https://byuzbasi.github.io/softridge/reference/hr_sure.md)
  : Stein unbiased risk estimate and degrees of freedom of the STR fit
- [`sure_hr()`](https://byuzbasi.github.io/softridge/reference/sure_hr.md)
  : STR estimator tuned by the Stein unbiased risk estimate

## Threshold and ridge primitives

- [`soft_threshold()`](https://byuzbasi.github.io/softridge/reference/soft_threshold.md)
  : Soft-thresholding operator
- [`ridge_path()`](https://byuzbasi.github.io/softridge/reference/ridge_path.md)
  : Ridge regression path computed through the singular value
  decomposition
