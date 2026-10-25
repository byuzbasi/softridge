#' Soft-thresholding operator
#'
#' Componentwise soft thresholding, \code{sign(z) * pmax(abs(z) - gamma, 0)}.
#'
#' @param z numeric vector or matrix.
#' @param gamma nonnegative threshold.
#' @return An object of the same shape as \code{z}.
#' @examples
#' soft_threshold(c(-2, -0.5, 0.3, 1.7), 1)
#' @export
soft_threshold <- function(z, gamma) {
  if (!length(gamma) || any(!is.finite(gamma)) || any(gamma < 0))
    stop("gamma must be finite and nonnegative")
  cpp_soft(z, gamma)
}

## Centre (and optionally scale) the design; returns the transformed matrix
## together with the centring and scaling constants so that coefficients can
## be mapped back to the original scale.
.prep_design <- function(X, y, intercept, standardize) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (!is.numeric(X) || any(!is.finite(X)) || any(!is.finite(y)))
    stop("X and y must be finite numeric data; no automatic imputation is performed")
  if (length(y) != n) stop("length(y) must equal nrow(X)")
  xm <- if (intercept) colMeans(X) else rep(0, ncol(X))
  Xc <- sweep(X, 2, xm, "-")
  xs <- if (standardize) sqrt(colSums(Xc^2) / n) else rep(1, ncol(X))
  if (any(xs == 0)) stop("design has constant (zero-variance) columns")
  Xs <- sweep(Xc, 2, xs, "/")
  ym <- if (intercept) mean(y) else 0
  list(X = Xs, y = y - ym, xm = xm, xs = xs, ym = ym)
}

## Thin singular value decomposition of the prepared design. Only singular
## values above a relative tolerance are kept, so that V spans the row space
## of X exactly (up to numerical precision).
.svd_design <- function(X, tol = 1e-10) {
  cpp_svd(X, tol)
}

#' Ridge regression path computed through the singular value decomposition
#'
#' Computes \eqn{(X^T X + k I_p)^{-1} X^T y = X^T (X X^T + k I_n)^{-1} y} for
#' every value in \code{k}. The singular value decomposition is computed once,
#' so each additional value of \code{k} costs \eqn{O(p r)} operations with
#' \eqn{r = rank(X)}.
#'
#' @param X prepared (centred/scaled) design matrix.
#' @param y prepared response.
#' @param k positive ridge parameters.
#' @param sv optional precomputed output of the internal SVD helper.
#' @return A \code{p x length(k)} matrix of ridge coefficients.
#' @keywords internal
ridge_path <- function(X, y, k, sv = NULL) {
  if (is.null(sv)) sv <- .svd_design(X)
  if (!length(k) || any(!is.finite(k)) || any(k <= 0)) stop("k must be positive and finite")
  cpp_ridge(sv$u, sv$d, sv$v, drop(y), k)
}

#' Soft-thresholded ridge (STR) estimator
#'
#' Fits the soft-thresholded ridge estimator
#' \eqn{\hat\theta_j = sign(\hat\beta_{k,j})(|\hat\beta_{k,j}| - \gamma)_+},
#' where \eqn{\hat\beta_k} is the ridge estimator with ridge parameter \code{k}.
#' The estimator has an explicit form and is valid for \eqn{p > n}.
#'
#' @param X \code{n x p} design matrix.
#' @param y response vector of length \code{n}.
#' @param k ridge parameter(s), positive.
#' @param gamma threshold(s), nonnegative. Thresholds apply on the scale of the
#'   standardized design when \code{standardize = TRUE}.
#' @param intercept logical, fit an unpenalized intercept by centring.
#' @param standardize logical, scale columns of \code{X} to unit variance.
#' @param scale For \code{hr}, logical training least-squares rescaling (SSTR-LS).
#'   For \code{hr_risk}, a fixed positive numeric scale; not a random fitted scale.
#' @param nonneg logical, use the nonnegative variant
#'   \eqn{(\hat\beta_{k,j} - \gamma)_+}.
#' @return An object of class \code{"hr"} containing the coefficient array
#'   \code{beta} of dimension \code{p x length(k) x length(gamma)} on the
#'   original scale, the intercepts \code{a0}, and the tuning grids.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(50 * 200), 50, 200)
#' y <- X[, 1:3] %*% c(3, -2, 1.5) + rnorm(50)
#' fit <- hr(X, y, k = 10, gamma = 0.1)
#' which(coef(fit)[-1] != 0)
#' @export
hr <- function(X, y, k, gamma, intercept = TRUE, standardize = TRUE,
               nonneg = FALSE, scale = FALSE) {
  if (length(scale)!=1L || is.na(scale) || !is.logical(scale)) stop("scale must be TRUE or FALSE")
  if (any(k <= 0)) stop("k must be positive")
  if (!length(gamma) || any(!is.finite(gamma)) || any(gamma < 0))
    stop("gamma must be finite and nonnegative")
  pr <- .prep_design(X, y, intercept, standardize)
  sv <- .svd_design(pr$X)
  B <- ridge_path(pr$X, pr$y, k, sv)
  path <- cpp_hr_path(B, gamma, pr$xs, pr$xm, pr$ym, nonneg, pr$X, drop(pr$y), scale)
  beta <- path$beta; a0 <- path$a0
  structure(list(beta = beta, a0 = a0, k = k, gamma = gamma,
                 ridge = B / pr$xs, nonneg = nonneg, scale = path$scale,
                 call = match.call()), class = "hr")
}

#' @param object a fitted \code{"hr"} object.
#' @param ik,ig indices of the ridge parameter and threshold to use.
#' @param newx matrix of new predictor values.
#' @param ... unused.
#' @rdname hr
#' @export
coef.hr <- function(object, ik = 1, ig = 1, ...) {
  c(object$a0[ik, ig], object$beta[, ik, ig])
}

#' @rdname hr
#' @export
predict.hr <- function(object, newx, ik = 1, ig = 1, ...) {
  drop(object$a0[ik, ig] + as.matrix(newx) %*% object$beta[, ik, ig])
}

#' Exact finite-sample moments and risk of the STR estimator
#'
#' Under \eqn{Y = X\beta + \epsilon}, \eqn{\epsilon \sim N_n(0, \sigma^2 I_n)}
#' and a fixed design, returns the exact mean, second moment and mean squared
#' error of each coordinate of the STR estimator as an estimator of the
#' projection \eqn{\theta = P_X \beta}, and the total risk
#' \eqn{E\|\hat\theta - \theta\|_2^2}. No centring or scaling is applied:
#' \code{X} is used as supplied.
#'
#' @param X design matrix.
#' @param beta true coefficient vector.
#' @param sigma error standard deviation.
#' @param k ridge parameter.
#' @param gamma threshold.
#' @param scale For \code{hr}, logical training least-squares rescaling (SSTR-LS).
#'   For \code{hr_risk}, a fixed positive numeric scale; not a random fitted scale.
#' @param nonneg logical, nonnegative variant.
#' @return A list with components \code{mean}, \code{second}, \code{mse},
#'   \code{risk}, \code{theta}, \code{m} (mean of the ridge estimator) and
#'   \code{s} (its coordinatewise standard deviations).
#' @importFrom stats dnorm pnorm
#' @export
hr_risk <- function(X, beta, sigma, k, gamma, nonneg = FALSE, scale = 1) {
  .validate_fixed_scale(scale)
  if(!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
     !is.numeric(beta) || length(beta)!=ncol(X) || any(!is.finite(beta)))stop("finite design and matching beta required")
  if(length(sigma)!=1L || !is.finite(sigma) || sigma<=0 ||
     length(k)!=1L || !is.finite(k) || k<=0 ||
     length(gamma)!=1L || !is.finite(gamma) || gamma<0)stop("sigma and k must be positive; gamma must be nonnegative scalars")
  cpp_risk(X, beta, sigma, k, gamma, nonneg, scale)
}

#' Stein unbiased risk estimate and degrees of freedom of the STR fit
#'
#' For a fixed design without intercept, returns the degrees of freedom
#' \eqn{\sum_{j \in \hat S} (A_k)_{jj}}, with
#' \eqn{A_k = (X^T X + k I)^{-1} X^T X}, and the unbiased estimate of the
#' prediction risk \eqn{E\|X\hat\theta - X\beta\|^2} for fixed tuning and known error variance.
#' Plug-in variance estimates or tuning by minimizing SURE do not retain that guarantee. When
#' \code{intercept = TRUE} the design and response are centred and one
#' degree of freedom is added for the intercept.
#'
#' @param X design matrix.
#' @param y response.
#' @param scale Fixed positive numeric scale of the slope fit. Intercept df remains one.
#'   These formulas do not apply to a scale estimated from the same response.
#' @param sigma2 error variance (known or estimated).
#' @param k ridge parameters.
#' @param gamma thresholds.
#' @param intercept logical.
#' @param standardize logical.
#' @return A list with matrices \code{sure} and \code{df} of dimension
#'   \code{length(k) x length(gamma)}.
#' @export
hr_sure <- function(X, y, sigma2, k, gamma, intercept = TRUE,
                    standardize = TRUE, scale = 1) {
  .validate_fixed_scale(scale)
  if(length(sigma2)!=1L || !is.finite(sigma2) || sigma2<=0)stop("sigma2 must be positive and finite")
  pr <- .prep_design(X, y, intercept, standardize)
  sv <- .svd_design(pr$X)
  B <- ridge_path(pr$X, pr$y, k, sv)
  if (!length(gamma) || any(!is.finite(gamma)) || any(gamma < 0))
    stop("gamma must be finite and nonnegative")
  G <- matrix(rep(gamma, each=length(k)), nrow=length(k))
  out <- cpp_sure_grid(pr$X, drop(pr$y), B, sv$v, sv$d, k, G, sigma2, intercept, scale)
  list(sure = out$sure, df = out$df, k = k, gamma = gamma, scale_min = out$scale_min)
}

.validate_fixed_scale <- function(scale) {
  if(length(scale)!=1L || !is.numeric(scale) || !is.finite(scale) || scale<=0)
    stop("fixed scale must be one positive finite number; data-dependent scale has no exact SURE guarantee")
}
.ls_scale <- function(X,y,Th) cpp_ls_scale(X,drop(y),as.matrix(Th))
