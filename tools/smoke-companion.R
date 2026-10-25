library(softridge)
library(scaleCorrection)
set.seed(25)
X <- matrix(rnorm(60 * 100), 60, 100)
y <- drop(X[, 1:3] %*% c(2, -1.5, 1) + rnorm(60))
grid <- hr_grid(X, y, nk = 3, ng = 5)
foldid <- rep(1:3, length.out = nrow(X))
pair <- corrected_pair(X, y,
  grid = list(k = grid$k, fraction = grid$fraction),
  foldid = foldid, unit_id = "showcase-25")
stopifnot(length(pair$SSTR_c0$coef) == 101L,
  length(pair$Ridge_c0$coef) == 101L,
  all(is.finite(pair$cv_error)), pair$Ridge_c0$gamma == 0,
  pair$SSTR_c0$scale >= 0, pair$Ridge_c0$scale >= 0)
print(pair$SSTR_c0[c("k", "gamma", "scale", "size")])
print(pair$Ridge_c0[c("k", "gamma", "scale", "size")])
cat("BOUNDED COMPANION SMOKE PASSED; one 60-by-100 example; no study rerun.\n")
