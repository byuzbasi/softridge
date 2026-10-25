#' Default tuning grids for the STR estimator
#'
#' The ridge grid is geometric between \code{kmin_ratio * d1^2} and
#' \code{d1^2}, where \code{d1} is the largest singular value of the prepared
#' design. For each ridge value the threshold grid is
#' \code{c * max|ridge|}, \code{c} in \code{[0, 0.99]}.
#'
#' @param X design matrix.
#' @param y response.
#' @param nk number of ridge values.
#' @param ng number of threshold values.
#' @param kmin_ratio smallest ridge value relative to the squared leading
#'   singular value.
#' @param intercept,standardize passed to the design preparation.
#' @return A list with the vector \code{k} and the matrix \code{gamma}
#'   (\code{nk x ng}), one row of thresholds per ridge value, and
#'   \code{fraction}, recomputed on the training response inside each CV fold.
#' @export
hr_grid <- function(X, y, nk = 12, ng = 40, kmin_ratio = 1e-4,
                    intercept = TRUE, standardize = TRUE) {
  if(any(!is.finite(c(nk,ng,kmin_ratio))) || length(nk)!=1L || length(ng)!=1L ||
     nk<1 || ng<1 || nk!=as.integer(nk) || ng!=as.integer(ng) ||
     length(kmin_ratio)!=1L || kmin_ratio<=0 || kmin_ratio>1)stop("invalid grid dimensions or ridge ratio")
  pr <- .prep_design(X, y, intercept, standardize)
  sv <- .svd_design(pr$X)
  k <- sv$d[1]^2 * exp(seq(log(kmin_ratio), 0, length.out = nk))
  B <- ridge_path(pr$X, pr$y, k, sv)
  cgrid <- seq(0, 0.99, length.out = ng)
  gamma <- t(vapply(seq_len(nk), function(i) cgrid * max(abs(B[, i])), numeric(ng)))
  dim(gamma) <- c(nk, ng)
  list(k = k, gamma = gamma, fraction = cgrid)
}

## Explicit, caller-supplied absolute grids remain supported. Default grids
## carry fractions, whose response-dependent scale is fitted inside each fold.
.fold_gamma <- function(grid, B, i) {
  if (is.null(grid$fraction)) grid$gamma[i, ] else
    grid$fraction * max(abs(B[, i]))
}

.validate_folds <- function(foldid, n) {
  if (length(foldid) != n || anyNA(foldid) ||
      any(foldid != as.integer(foldid)) || min(foldid) != 1L ||
      max(foldid) < 2L || !identical(sort(unique(as.integer(foldid))),
                                  seq_len(max(foldid))))
    stop("foldid must contain consecutive integers starting at 1 and at least two folds")
  if (any(n - tabulate(foldid) < 3L)) stop("each training fold needs at least 3 observations")
  invisible(TRUE)
}

## Moore-Penrose least squares on the prepared design, using the same
## relative singular-value tolerance as the ridge decomposition.
.mp_refit <- function(X, y) {
  cpp_mp_refit(X, drop(y))
}

#' Cross-validated STR estimator
#'
#' Selects the ridge parameter and the threshold by K-fold cross-validation of
#' the prediction error. With \code{refit = TRUE}, the variables retained by
#' thresholding are refitted by least squares and the pair (ridge parameter,
#' model size) is selected instead; model size \code{m} corresponds to the
#' data-dependent threshold equal to the \code{(m+1)}-th largest absolute
#' ridge coefficient when magnitudes are distinct. Ties use column order.
#' Refit size is capped by the smallest training-fold rank, and refits use
#' the Moore-Penrose solution on the prepared design (relative tolerance 1e-10).
#'
#' @param X design matrix.
#' @param y response.
#' @param nfolds number of folds.
#' @param grid output of \code{hr_grid}; computed on the full data if missing.
#' @param scale Logical, estimate a least-squares scale inside each training fold
#'   and on the final training sample (SSTR-LS). Incompatible with refit.
#' @param refit logical, least-squares refit on the selected variables.
#' @param max_size largest model size considered when \code{refit = TRUE}.
#' @param foldid optional fold assignment.
#' @param intercept,standardize passed to the design preparation.
#' @return A list with the selected tuning values, the cross-validation error
#'   surface, the final coefficient vector \code{coef} (intercept first) and
#'   the selected index set \code{selected}.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 150), 60, 150)
#' y <- X[, 1:3] %*% c(3, -2, 1.5) + rnorm(60)
#' cv <- cv_hr(X, y, nfolds = 5)
#' cv$selected
#' @name cv_hr
NULL

.cv_threshold <- function(X, y, nfolds = 10, grid = NULL, refit = FALSE,
                  max_size = NULL, foldid = NULL, intercept = TRUE,
                  standardize = TRUE, scale = FALSE, hard = FALSE) {
  if(length(scale)!=1L || is.na(scale) || !is.logical(scale)) stop("scale must be TRUE or FALSE")
  if(refit && (scale || hard)) stop("use either refit = TRUE or scale = TRUE")
  X <- as.matrix(X); y <- drop(y); n <- nrow(X); p <- ncol(X)
  if (is.null(grid)) grid <- hr_grid(X, y, intercept = intercept,
                                     standardize = standardize)
  if (is.null(foldid)) foldid <- sample(rep(seq_len(nfolds), length.out = n))
  .validate_folds(foldid, n)
  nfolds <- max(foldid)
  k <- grid$k
  if (refit) {
    if (is.null(max_size)) max_size <- min(floor(0.5 * n), p)
    if (length(max_size) != 1L || !is.finite(max_size) || max_size < 1 ||
        max_size != as.integer(max_size)) stop("max_size must be a positive integer")
    # A common candidate set across folds, bounded by their design ranks.
    fold_ranks <- vapply(seq_len(nfolds), function(f) {
      tr <- foldid != f
      pr <- .prep_design(X[tr, , drop = FALSE], y[tr], intercept, standardize)
      length(.svd_design(pr$X)$d)
    }, integer(1))
    requested_max_size <- max_size
    max_size <- min(max_size, p, fold_ranks)
    if (max_size < 1L) stop("no estimable refit candidates")
    err <- matrix(0, length(k), max_size)
  } else {
    err <- matrix(0, length(k), ncol(grid$gamma))
  }
  for (f in seq_len(nfolds)) {
    tr <- foldid != f
    pr <- .prep_design(X[tr, , drop = FALSE], y[tr], intercept, standardize)
    sv <- .svd_design(pr$X)
    B <- ridge_path(pr$X, pr$y, k, sv)
    Xte <- sweep(sweep(X[!tr, , drop = FALSE], 2, pr$xm, "-"), 2, pr$xs, "/")
    yte <- y[!tr] - pr$ym
    G <- if (refit) matrix(0, length(k), 0) else
      t(vapply(seq_along(k), function(i) .fold_gamma(grid, B, i), numeric(ncol(grid$gamma))))
    if (!refit) dim(G) <- c(length(k), ncol(grid$gamma))
    err <- err + cpp_cv_scores(pr$X, drop(pr$y), Xte, drop(yte), B, G,
                                if (refit) as.integer(max_size) else 0L, scale, hard)
  }
  err <- err / n
  best <- which(err == min(err), arr.ind = TRUE)[1, ]
  pr <- .prep_design(X, y, intercept, standardize)
  sv <- .svd_design(pr$X)
  b <- drop(ridge_path(pr$X, pr$y, k[best[1]], sv))
  cc <- 1
  if (refit) {
    S <- sort(order(abs(b), decreasing = TRUE)[seq_len(best[2])])
    th <- numeric(p)
    refit_info <- .mp_refit(pr$X[, S, drop = FALSE], pr$y)
    th[S] <- refit_info$coef
    gamma_sel <- if (best[2] < p) sort(abs(b), decreasing = TRUE)[best[2] + 1] else 0
  } else {
    gamma_sel <- if (is.null(grid$fraction)) grid$gamma[best[1], best[2]] else
      grid$fraction[best[2]] * max(abs(b))
    th <- if(hard) b*(abs(b)>gamma_sel) else soft_threshold(b, gamma_sel)
    if(scale) { cc <- .ls_scale(pr$X,pr$y,th)[1]; th <- cc*th }
  }
  beta <- th / pr$xs
  a0 <- pr$ym - sum(beta * pr$xm)
  list(k = k[best[1]], gamma = gamma_sel, scale = cc, threshold_rule = if(hard) "hard" else "soft",
       size = if (refit) best[2] else sum(th != 0),
       cv_error = err, min_error = min(err), grid = grid, foldid = foldid,
       coef = c(a0, beta), selected = which(beta != 0), refit = refit,
       refit_diagnostics = if (refit) list(rank = refit_info$rank,
         requested_max_size = requested_max_size, effective_max_size = max_size,
         fold_ranks = fold_ranks, tie_rule = "original column order") else NULL)
}

#' STR estimator tuned by the Stein unbiased risk estimate
#'
#' Minimizes the plug-in SURE criterion over the grid. Its selected minimum
#' is not an unbiased estimate of the risk of the tuned estimator.
#'
#' @param X design matrix.
#' @param y response.
#' @param sigma2 error variance or an estimate of it.
#' @param grid output of \code{hr_grid}; computed if missing.
#' @param intercept,standardize passed to the design preparation.
#' @return A list with the selected tuning values, the fitted coefficients
#'   (intercept first), the selected set and the degrees of freedom.
#' @export
sure_hr <- function(X, y, sigma2, grid = NULL, intercept = TRUE,
                    standardize = TRUE) {
  if(length(sigma2)!=1L || !is.finite(sigma2) || sigma2<=0)stop("sigma2 must be positive and finite")
  X <- as.matrix(X); y <- drop(y)
  if (is.null(grid)) grid <- hr_grid(X, y, intercept = intercept,
                                     standardize = standardize)
  pr <- .prep_design(X, y, intercept, standardize)
  sv <- .svd_design(pr$X)
  B <- ridge_path(pr$X, pr$y, grid$k, sv)
  surface <- cpp_sure_grid(pr$X, drop(pr$y), B, sv$v, sv$d, grid$k,
                           grid$gamma, sigma2, intercept)
  best <- c(Inf, NA, NA); dfb <- NA
  # Preserve the original k-first then gamma tie-breaking order.
  for (i in seq_along(grid$k)) {
    j <- which.min(surface$sure[i, ])
    if (surface$sure[i, j] < best[1]) {
      best <- c(surface$sure[i, j], i, j); dfb <- surface$df[i, j]
    }
  }
  kk <- grid$k[best[2]]; gg <- grid$gamma[best[2], best[3]]
  b <- drop(ridge_path(pr$X, pr$y, kk, sv))
  beta <- soft_threshold(b, gg) / pr$xs
  a0 <- pr$ym - sum(beta * pr$xm)
  list(k = kk, gamma = gg, sure = best[1], df = dfb, coef = c(a0, beta),
       selected = which(beta != 0),
       risk_status = "plug-in tuning criterion; not unbiased risk of the tuned fit")
}

#' @rdname cv_hr
#' @export
cv_hr <- function(X,y,nfolds=10,grid=NULL,refit=FALSE,max_size=NULL,foldid=NULL,
                  intercept=TRUE,standardize=TRUE,scale=FALSE) {
  .cv_threshold(X,y,nfolds,grid,refit,max_size,foldid,intercept,standardize,scale,FALSE)
}
.cv_hard_ridge <- function(X,y,grid,foldid,scale=FALSE) {
  .cv_threshold(X,y,grid=grid,foldid=foldid,scale=scale,hard=TRUE)
}
