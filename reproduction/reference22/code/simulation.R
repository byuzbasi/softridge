###############################################################################
## simulation.R
##
## Simulation study for "High-dimensional ridge shrinkage and selection by
## soft thresholding" (Arashi and Yuzbasi).
##
## Five parts:
##   Study 1  Exact finite-sample risk of the HR estimator (Theorem 1) and the
##            unbiasedness of the Stein risk estimate (Theorem 2), checked
##            against Monte Carlo averages for a fixed p > n design, and the
##            exact risk-minimizing scale factor (Proposition 2).
##   Study 2  Comparison with ridge, lasso, elastic net, adaptive lasso, SCAD,
##            MCP, broken adaptive ridge and the hard-thresholded ridge of
##            Shao and Deng (2012), all tuned by 10-fold cross-validation on
##            common folds, for a sparse regression vector and weakly
##            correlated predictors.
##   Study 3  Computing time as a function of (n, p).
##   Study 4  The same comparison for block-correlated predictors, with a
##            regression vector that is constant within blocks (scenarios 1 to
##            4) and four departures from that setting (scenarios 5 to 8).
##   Study 5  Low signal-to-noise ratio, SNR = 0.2, for the designs of
##            Studies 2 and 4, where the least-squares scale overfits.
##
## Data source: all data are generated inside this script; no external data.
## Requirements: R (>= 4.1), glmnet, ncvreg, and the accompanying package
##   softridge (install with R CMD INSTALL softridge_0.4.0.tar.gz).
## Output: results/*.rds, results/*.csv and figures/*.pdf, relative to the
##   working directory.
###############################################################################

suppressPackageStartupMessages({
  library(softridge)   # HR estimator, exact risk, SURE, cross-validation
  library(glmnet)      # ridge, lasso, elastic net, adaptive lasso
  library(ncvreg)      # SCAD and MCP
  library(parallel)    # replications in parallel
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
NCORES <- max(1, min(2, detectCores()))
RNGkind("L'Ecuyer-CMRG")
set.seed(20260930)

## Colours: one highlight hue for the proposed method, grey for competitors.
COL_HR   <- "#2a78d6"
COL_ALT  <- "#eb6834"
COL_GREY <- "#8c8c8c"

###############################################################################
## Helpers
###############################################################################

## Covariance matrices: first-order autoregressive and compound symmetry.
ar1 <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))
cs  <- function(p, rho) { S <- matrix(rho, p, p); diag(S) <- 1; S }

## Draw an n x p Gaussian design with covariance Sigma, using a Cholesky
## factor computed once per scenario.
rdesign <- function(n, R) matrix(rnorm(n * nrow(R)), n) %*% R

## Least-squares scale of each column t of Th: the value c >= 0 minimizing
## ||y - c X t||^2 (one if X t = 0).
ls_scale <- function(X, y, Th) {
  Fit <- X %*% as.matrix(Th); den <- colSums(Fit^2)
  ifelse(den > 0, pmax(drop(crossprod(Fit, y)) / pmax(den, .Machine$double.xmin), 0), 1)
}

## Optimism-corrected scale (Proposition 2(ii)): the least-squares scale
## reduced by sigma2 * df / ||X t||^2, where df is the sum of the ridge
## leverages ad = diag(A_k) over the nonzero coordinates of t; truncated at zero.
corrected_scale <- function(X, y, Th, ad, sigma2) {
  Th <- as.matrix(Th); Fit <- X %*% Th; den <- colSums(Fit^2)
  num <- drop(crossprod(Fit, y)) - sigma2 * colSums((Th != 0) * ad)
  ifelse(den > 0, pmax(num / pmax(den, .Machine$double.xmin), 0), 1)
}

## Hard-thresholded ridge estimator of Shao and Deng (2012), tuned by
## K-fold cross-validation over the same (k, gamma) grid as HR, so that the
## only difference between the two is the thresholding rule. With
## scale = TRUE the hard-thresholded vector is multiplied by the corrected
## scale, exactly as for the scaled HR estimator; sigma2 is one value or the
## fold-specific vector of fold_variance(), used in the same way as in cv_hr().
cv_hard_ridge <- function(X, y, grid, foldid, scale = FALSE, sigma2 = NULL) {
  n <- nrow(X); err <- matrix(0, length(grid$k), ncol(grid$gamma))
  K <- max(foldid); if (scale && length(sigma2) == 1) sigma2 <- rep(sigma2, K + 1)
  prep <- function(Xt, yt) {
    xm <- colMeans(Xt); Xc <- sweep(Xt, 2, xm)
    xs <- sqrt(colSums(Xc^2) / nrow(Xt))
    list(X = sweep(Xc, 2, xs, "/"), y = yt - mean(yt), xm = xm, xs = xs, ym = mean(yt))
  }
  ridge <- function(Xs, ys, k) {
    s <- svd(Xs); keep <- s$d > 1e-10 * s$d[1]
    u <- s$u[, keep]; d <- s$d[keep]; v <- s$v[, keep]; uty <- crossprod(u, ys)
    B <- sapply(k, function(kk) drop(v %*% (d / (d^2 + kk) * uty)))
    attr(B, "ad") <- sapply(k, function(kk) drop((v^2) %*% (d^2 / (d^2 + kk))))
    B
  }
  for (f in seq_len(max(foldid))) {
    tr <- foldid != f; pr <- prep(X[tr, ], y[tr]); B <- ridge(pr$X, pr$y, grid$k)
    Xte <- sweep(sweep(X[!tr, , drop = FALSE], 2, pr$xm), 2, pr$xs, "/")
    yte <- y[!tr] - pr$ym
    for (i in seq_along(grid$k)) {
      Th <- sapply(grid$gamma[i, ], function(g) B[, i] * (abs(B[, i]) > g))
      if (scale) Th <- sweep(Th, 2, corrected_scale(pr$X, pr$y, Th, attr(B, "ad")[, i], sigma2[f]), "*")
      err[i, ] <- err[i, ] + colSums((yte - Xte %*% Th)^2)
    }
  }
  best <- which(err == min(err), arr.ind = TRUE)[1, ]
  pr <- prep(X, y); B <- ridge(pr$X, pr$y, grid$k[best[1]]); b <- drop(B)
  th <- b * (abs(b) > grid$gamma[best[1], best[2]])
  if (scale) th <- corrected_scale(pr$X, pr$y, th, drop(attr(B, "ad")), sigma2[K + 1]) * th
  th <- th / pr$xs
  c(pr$ym - sum(th * pr$xm), th)
}

## Broken adaptive ridge (Dai et al., 2018): iteratively reweighted ridge,
## beta^(t+1) = (X'X + lambda D_t)^{-1} X'y with D_t = diag(1/beta_j^(t)^2),
## started from the ridge fit with parameter xi. Computed in the
## n-dimensional form W X'(X W X' + lambda I_n)^{-1} y with W = diag(beta^(t)^2),
## which is exact and keeps coefficients that have reached zero at zero.
## lambda is chosen by 10-fold cross-validation on the common folds.
bar_path <- function(Xs, yc, lambdas, xi = 1, maxit = 100, tol = 1e-6) {
  n <- nrow(Xs)
  b0 <- drop(crossprod(Xs, solve(tcrossprod(Xs) + xi * diag(n), yc)))
  sapply(lambdas, function(lam) {
    b <- b0
    for (it in seq_len(maxit)) {
      act <- which(abs(b) > 1e-8)
      if (length(act) == 0) break
      Xa <- Xs[, act, drop = FALSE]; w <- b[act]^2
      bn <- numeric(length(b))
      bn[act] <- w * drop(crossprod(Xa, solve(Xa %*% (w * t(Xa)) + lam * diag(n), yc)))
      conv <- max(abs(bn - b)) < tol; b <- bn
      if (conv) break
    }
    b[abs(b) < 1e-6] <- 0
    b
  })
}
cv_bar <- function(X, y, foldid, lambdas = exp(seq(log(0.01), log(100), length.out = 15))) {
  prep <- function(Xt, yt) {
    xm <- colMeans(Xt); Xc <- sweep(Xt, 2, xm)
    xs <- sqrt(colSums(Xc^2) / nrow(Xt))
    list(X = sweep(Xc, 2, xs, "/"), y = yt - mean(yt), xm = xm, xs = xs, ym = mean(yt))
  }
  err <- numeric(length(lambdas))
  for (f in seq_len(max(foldid))) {
    tr <- foldid != f; pr <- prep(X[tr, ], y[tr])
    B <- bar_path(pr$X, pr$y, lambdas)
    Xte <- sweep(sweep(X[!tr, , drop = FALSE], 2, pr$xm), 2, pr$xs, "/")
    err <- err + colSums((y[!tr] - pr$ym - Xte %*% B)^2)
  }
  pr <- prep(X, y)
  b <- drop(bar_path(pr$X, pr$y, lambdas[which.min(err)])) / pr$xs
  c(pr$ym - sum(b * pr$xm), b)
}

## Error variance estimate from the cross-validated lasso (Reid et al., 2016):
## residual sum of squares divided by n minus the size of the lasso model.
sigma2_lasso <- function(fit, X, y) {
  b <- as.numeric(coef(fit, s = "lambda.min"))
  s <- sum(b[-1] != 0)
  sum((y - b[1] - X %*% b[-1])^2) / max(1, length(y) - s - 1)
}

## The same estimate computed from a given sample alone, with the lasso tuned
## by cross-validation over the supplied folds. fold_variance() applies it to
## the training part of every fold, with the remaining folds as inner folds,
## so that the held-out responses never enter the scale used to predict them.
sigma2_cvlasso <- function(X, y, foldid)
  sigma2_lasso(cv.glmnet(X, y, alpha = 1, foldid = foldid), X, y)

## Fit every method on one data set; returns a list of coefficient vectors
## (intercept first) and fitting times in seconds.
fit_all <- function(X, y, nfolds = 10) {
  n <- nrow(X); p <- ncol(X)
  foldid <- sample(rep(seq_len(nfolds), length.out = n))
  out <- list(); tm <- c()
  tic <- function() proc.time()[3]

  ## Ridge regression over a long penalty path. The default path of glmnet
  ## spans two orders of magnitude when p > n and can end before the minimum
  ## of the cross-validation curve; six orders of magnitude are used here.
  t0 <- tic(); f <- cv.glmnet(X, y, alpha = 0, foldid = foldid,
                              lambda.min.ratio = 1e-6, nlambda = 200)
  out$Ridge <- as.numeric(coef(f, s = "lambda.min")); tm["Ridge"] <- tic() - t0

  t0 <- tic(); flas <- cv.glmnet(X, y, alpha = 1, foldid = foldid)
  out$Lasso <- as.numeric(coef(flas, s = "lambda.min")); tm["Lasso"] <- tic() - t0

  ## Elastic net: the mixing parameter is chosen on the same folds, from a
  ## set that includes nearly ridge-like mixtures.
  t0 <- tic(); cvs <- lapply(c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8), function(a)
    cv.glmnet(X, y, alpha = a, foldid = foldid))
  j <- which.min(sapply(cvs, function(f) min(f$cvm)))
  out$ENet <- as.numeric(coef(cvs[[j]], s = "lambda.min")); tm["ENet"] <- tic() - t0

  ## Adaptive lasso with weights 1/|lasso|; variables dropped by the initial
  ## lasso receive an infinite weight and are excluded.
  t0 <- tic(); b0 <- out$Lasso[-1]; keep <- which(b0 != 0)
  if (length(keep) >= 2) {
    fa <- cv.glmnet(X[, keep, drop = FALSE], y, alpha = 1, foldid = foldid,
                    penalty.factor = 1 / abs(b0[keep]))
    ca <- as.numeric(coef(fa, s = "lambda.min"))
    out$ALasso <- c(ca[1], replace(numeric(p), keep, ca[-1]))
  } else out$ALasso <- out$Lasso
  tm["ALasso"] <- tic() - t0 + tm["Lasso"]

  t0 <- tic(); fs <- cv.ncvreg(X, y, penalty = "SCAD", fold = foldid)
  out$SCAD <- as.numeric(coef(fs)); tm["SCAD"] <- tic() - t0
  t0 <- tic(); fm <- cv.ncvreg(X, y, penalty = "MCP", fold = foldid)
  out$MCP <- as.numeric(coef(fm)); tm["MCP"] <- tic() - t0

  t0 <- tic(); out$BAR <- cv_bar(X, y, foldid)
  tm["BAR"] <- tic() - t0

  grid <- hr_grid(X, y)
  t0 <- tic(); out$HardRidge <- cv_hard_ridge(X, y, grid, foldid); tm["HardRidge"] <- tic() - t0
  t0 <- tic(); cvh <- cv_hr(X, y, grid = grid, foldid = foldid)
  out$HR <- cvh$coef; tm["HR"] <- tic() - t0
  t0 <- tic(); s2 <- sigma2_lasso(flas, X, y)
  out$HR_SURE <- sure_hr(X, y, s2, grid = grid)$coef; tm["HR_SURE"] <- tic() - t0
  t0 <- tic(); out$HR_OLS <- cv_hr(X, y, grid = grid, foldid = foldid, refit = TRUE,
                                   max_size = 40)$coef
  tm["HR_OLS"] <- tic() - t0
  t0 <- tic(); s2f <- fold_variance(X, y, foldid, sigma2_cvlasso)
  ## The full-data element repeats the lasso fit above on the same folds.
  stopifnot(abs(s2f[max(foldid) + 1] - s2) <= 1e-10 * s2)
  nw <- fit_scaled(X, y, foldid, grid, s2f); nw$time <- nw$time + (tic() - t0) * c(1, 1, 0, 1)
  list(coef = c(out, nw$coef), time = c(tm, nw$time))
}

## The scaled estimators. HR_S (rule T4): the soft-thresholded vector times
## the corrected scale, recomputed on every training set with the error
## variance estimated on that training set by the rule of (T2); s2 holds these
## fold-specific estimates followed by the full-data one. Ridge_S: the same with the threshold
## fixed at zero, a control that isolates the contribution of thresholding.
## HR_LS: the least-squares scale instead of the corrected one.
## HardRidge_S: the corrected scale applied to the hard-thresholded ridge.
fit_scaled <- function(X, y, foldid, grid, s2) {
  out <- list(); tm <- c(); tic <- function() proc.time()[3]
  t0 <- tic(); out$Ridge_S <- cv_hr(X, y, grid = grid, foldid = foldid, scale = TRUE,
                                    sigma2 = s2, zero_threshold = TRUE)$coef
  tm["Ridge_S"] <- tic() - t0
  t0 <- tic(); out$HR_S <- cv_hr(X, y, grid = grid, foldid = foldid, scale = TRUE,
                                 sigma2 = s2)$coef
  tm["HR_S"] <- tic() - t0
  t0 <- tic(); out$HR_LS <- cv_hr(X, y, grid = grid, foldid = foldid, scale = TRUE,
                                  scale_type = "ls")$coef
  tm["HR_LS"] <- tic() - t0
  t0 <- tic(); out$HardRidge_S <- cv_hard_ridge(X, y, grid, foldid, scale = TRUE, sigma2 = s2)
  tm["HardRidge_S"] <- tic() - t0
  list(coef = out, time = tm)
}

## Performance measures against the true beta for a design with mean zero and
## covariance Sigma: squared L2 estimation error, prediction error
## E(x'beta_hat + a0 - x'beta)^2 = a0^2 + (beta_hat - beta)' Sigma (beta_hat - beta),
## true and false positive rates, and model size.
evaluate <- function(cf, beta, Sigma) {
  a0 <- cf[1]; b <- cf[-1]; d <- b - beta; S <- beta != 0
  c(L2 = sqrt(sum(d^2)),
    PE = a0^2 + drop(crossprod(d, Sigma %*% d)),
    TPR = mean(b[S] != 0), FPR = mean(b[!S] != 0), Size = sum(b != 0))
}

###############################################################################
## Study 1: exact risk and SURE versus Monte Carlo
###############################################################################

study1 <- function(R = 20000) {
  set.seed(101)
  ## Two fixed designs with n = 50: p = 60 (projection close to the identity,
  ## so theta is nearly sparse) and p = 200 (theta dense).
  n <- 50; sigma <- 1
  kgrid <- c(1, 10, 100)
  ggrid <- seq(0, 0.3, by = 0.01)
  curves <- NULL
  for (p in c(60, 200)) {
    X <- rdesign(n, chol(ar1(p, 0.5)))
    X <- scale(X, center = TRUE, scale = FALSE)
    beta <- numeric(p); beta[c(1, 5, 10, 30, 50)] <- c(3, -2, 1.5, 1, -1)
    sv <- svd(X)
    ## Size and sparsity of the projection target theta = P beta (no random
    ## numbers are used here, so the simulation stream is unaffected).
    keep <- sv$d > 1e-10 * sv$d[1]
    theta <- drop(sv$v[, keep] %*% crossprod(sv$v[, keep], beta))
    cat(sprintf("p = %d: rank %d, ||theta - beta||/||beta|| = %.3f, #|theta_j| > 0.1 = %d, theta at signals: %s\n",
                p, sum(keep), sqrt(sum((theta - beta)^2) / sum(beta^2)), sum(abs(theta) > 0.1),
                paste(round(theta[beta != 0], 2), collapse = " ")))
    for (k in kgrid) {
      Hk <- sv$v %*% (sv$d / (sv$d^2 + k) * t(sv$u))          # p x n ridge map
      ad <- colSums(t(sv$v^2) * (sv$d^2 / (sv$d^2 + k)))       # diagonal of A_k
      E <- matrix(rnorm(n * R, sd = sigma), n, R)
      Y <- drop(X %*% beta) + E; B <- Hk %*% Y
      ex0 <- hr_risk(X, beta, sigma, k, 0)
      ## Right derivative of the exact risk at gamma = 0 (Proposition 1).
      b <- ex0$m - ex0$theta
      Delta <- sum(b * (2 * pnorm(ex0$m / ex0$s) - 1) + 2 * ex0$s * dnorm(ex0$m / ex0$s))
      for (g in ggrid) {
        ex <- hr_risk(X, beta, sigma, k, g)
        Th <- soft_threshold(B, g)
        loss <- colSums((Th - ex$theta)^2)
        fit <- X %*% Th
        ploss <- colSums((fit - drop(X %*% beta))^2)
        sure <- colSums((Y - fit)^2) - n * sigma^2 + 2 * sigma^2 * colSums((Th != 0) * ad)
        ## Risk at the risk-minimizing fixed scale (Proposition 2), exact and
        ## by Monte Carlo on the same replications.
        loss_c <- colSums((ex$scale_opt * Th - ex$theta)^2)
        curves <- rbind(curves, data.frame(p = p, k = k, gamma = g, exact = ex$risk,
                                           mc = mean(loss), mc_se = sd(loss) / sqrt(R),
                                           pred = mean(ploss), sure = mean(sure),
                                           sure_se = sd(sure - ploss) / sqrt(R),
                                           Delta = Delta, scale_opt = ex$scale_opt,
                                           exact_c = ex$risk_opt, mc_c = mean(loss_c),
                                           mc_c_se = sd(loss_c) / sqrt(R)))
      }
    }
  }
  saveRDS(curves, "results/study1.rds")
  write.csv(curves, "results/study1.csv", row.names = FALSE)

  pdf("figures/fig1_exact_risk.pdf", width = 7, height = 6.2)
  par(mfrow = c(2, 2), mar = c(4, 4.2, 1.8, 0.5), mgp = c(2.5, 0.7, 0), las = 1,
      cex.axis = 0.8, cex.lab = 0.9, cex.main = 0.9)
  cols <- c(COL_HR, COL_ALT, "#1baf7a"); pchs <- c(16, 17, 15)
  panel <- function(d, yv, pv, ylab, main) {
    plot(NA, xlim = range(ggrid), ylim = range(d[, c(yv, pv)]), xlab = expression(gamma),
         ylab = ylab, bty = "l", main = main)
    for (i in seq_along(kgrid)) {
      e <- d[d$k == kgrid[i], ]
      lines(e$gamma, e[[yv]], col = cols[i], lwd = 2)
      idx <- seq(1, nrow(e), by = 3)
      points(e$gamma[idx], e[[pv]][idx], col = cols[i], pch = pchs[i], cex = 0.75)
    }
    legend("topleft", legend = paste("k =", kgrid), col = cols, lwd = 2, pch = pchs,
           bty = "n", cex = 0.75)
  }
  panel(curves[curves$p == 60, ], "exact", "mc", "Risk for the projection", "A   p = 60")
  panel(curves[curves$p == 200, ], "exact", "mc", "Risk for the projection", "B   p = 200")
  panel(curves[curves$p == 60, ], "pred", "sure", "Prediction risk", "C   p = 60")
  panel(curves[curves$p == 200, ], "pred", "sure", "Prediction risk", "D   p = 200")
  dev.off()
  curves
}

## Summary of Study 1: agreement of exact and Monte Carlo risks, unbiasedness
## of SURE, and the risk-minimizing threshold for each (p, k).
summarise1 <- function(curves) {
  curves$z_risk <- (curves$mc - curves$exact) / curves$mc_se
  curves$z_sure <- (curves$sure - curves$pred) / curves$sure_se
  curves$z_scale <- (curves$mc_c - curves$exact_c) / curves$mc_c_se
  tab <- do.call(rbind, lapply(split(curves, list(curves$p, curves$k)), function(e) {
    i <- which.min(e$exact); j <- which.min(e$exact_c)
    data.frame(p = e$p[1], k = e$k[1], Delta = e$Delta[1], risk_ridge = e$exact[1],
               risk_min = e$exact[i], gamma_min = e$gamma[i],
               gain_pct = 100 * (1 - e$exact[i] / e$exact[1]),
               scale_ridge = e$scale_opt[1], risk_ridge_scaled = e$exact_c[1],
               risk_min_scaled = e$exact_c[j], gamma_min_scaled = e$gamma[j],
               scale_min_scaled = e$scale_opt[j],
               gain_scaled_pct = 100 * (1 - e$exact_c[j] / e$exact[1]),
               max_abs_z_scale = max(abs(e$z_scale)),
               pred_ridge = e$pred[1], pred_min = min(e$pred),
               gamma_pred_min = e$gamma[which.min(e$pred)],
               max_abs_z_risk = max(abs(e$z_risk)), max_abs_z_sure = max(abs(e$z_sure)))
  }))
  write.csv(tab, "results/study1_summary.csv", row.names = FALSE)
  tab
}

###############################################################################
## Study 2: comparison with competing methods
###############################################################################

study2 <- function(nrep = 100, n = 100) {
  scen <- expand.grid(p = c(200, 1000), cov = c("AR", "CS"), snr = c(1, 5),
                      stringsAsFactors = FALSE)
  res <- list()
  for (s in seq_len(nrow(scen))) {
    ## Each scenario has its own seed and its own results file, so that the
    ## study can be resumed after an interruption with identical results.
    fn <- sprintf("results/study2_scenario%d.rds", s)
    if (file.exists(fn)) { res[[s]] <- readRDS(fn); next }
    set.seed(2000 + s)
    p <- scen$p[s]
    Sigma <- if (scen$cov[s] == "AR") ar1(p, 0.5) else cs(p, 0.3)
    R <- chol(Sigma)
    beta <- numeric(p)
    beta[c(1, 2, 5, 10, 20, 50, 100, 150)] <- c(3, 1.5, 2, -2.5, 1, -1, 2, -1.5)
    sigma <- sqrt(drop(crossprod(beta, Sigma %*% beta)) / scen$snr[s])
    one <- function(r) {
      X <- rdesign(n, R); y <- drop(X %*% beta) + rnorm(n, sd = sigma)
      f <- fit_all(X, y)
      ev <- t(sapply(f$coef, evaluate, beta = beta, Sigma = Sigma))
      data.frame(scenario = s, p = p, cov = scen$cov[s], snr = scen$snr[s],
                 sigma = sigma, rep = r, method = rownames(ev), ev,
                 time = f$time[rownames(ev)], row.names = NULL)
    }
    t0 <- proc.time()[3]
    res[[s]] <- do.call(rbind, mclapply(seq_len(nrep), one, mc.cores = NCORES,
                                        mc.set.seed = TRUE))
    saveRDS(res[[s]], fn)
    cat(sprintf("scenario %d/%d (p=%d, %s, SNR=%g) done in %.0f s\n", s,
                nrow(scen), p, scen$cov[s], scen$snr[s], proc.time()[3] - t0))
  }
  out <- do.call(rbind, res)
  saveRDS(out, "results/study2.rds")
  out
}

## Methods reported in the tables and figures of the paper, in display order.
## The least-squares-scaled HR (HR_LS) and the scaled hard-thresholded ridge
## (HardRidge_S) are reported in the Supplementary Material.
METH <- c("Ridge", "Lasso", "ENet", "ALasso", "SCAD", "MCP", "BAR",
          "HardRidge", "HR", "HR_SURE", "HR_OLS", "Ridge_S", "HR_S")
LAB <- c("Ridge", "Lasso", "ENet", "ALasso", "SCAD", "MCP", "BAR",
         "TRidge", "HR", "HR-SURE", "HR-OLS", "Ridge-S", "HR-S")
## Additional rows of the supplementary tables.
METH_X <- c(METH, "HR_LS", "HardRidge_S")
LAB_X <- c(LAB, "HR-LS", "TRidge-S")

## Summary table: mean (standard error) over replications.
summarise2 <- function(out, file = "results/study2_summary.csv") {
  agg <- function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
  keys <- unique(out[, c("scenario", "p", "cov", "snr", "method")])
  tab <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- merge(out, keys[i, ])
    data.frame(keys[i, ], t(unlist(lapply(d[, c("L2", "PE", "TPR", "FPR", "Size", "time")], agg))))
  }))
  write.csv(tab, file, row.names = FALSE)
  tab
}

## LaTeX table for the paper: mean (standard error) of prediction error,
## L2 error, TPR, FPR and model size for every method and scenario.
table2 <- function(out, file = "results/table_study2.tex", meth = METH, lab = LAB,
                   header = function(d) sprintf("$p=%d$, %s, SNR $=%g$", d$p[1],
                     ifelse(d$cov[1] == "AR", "autoregressive", "compound symmetry"), d$snr[1])) {
  f <- function(v, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(v), sd(v) / sqrt(length(v)))
  lines <- c()
  for (s in sort(unique(out$scenario))) {
    d <- out[out$scenario == s, ]
    lines <- c(lines, sprintf("\\multicolumn{6}{l}{%s}\\\\", header(d)))
    for (i in seq_along(meth)) {
      e <- d[d$method == meth[i], ]
      lines <- c(lines, sprintf("%s & %s & %s & %.2f & %.3f & %.1f\\\\", lab[i], f(e$PE),
                                f(e$L2), mean(e$TPR), mean(e$FPR), mean(e$Size)))
    }
    lines <- c(lines, "\\addlinespace")
  }
  writeLines(lines, file)
  invisible(lines)
}
## Two tables for the paper, one per signal-to-noise ratio.
tables2 <- function(out) {
  table2(out[out$snr == 5, ], "results/table_study2_snr5.tex")
  table2(out[out$snr == 1, ], "results/table_study2_snr1.tex")
}

plot2 <- function(out, meth = METH, lab = LAB) {
  cols <- ifelse(grepl("^HR", meth), COL_HR, COL_GREY)
  pdf("figures/fig2_prediction_error.pdf", width = 7.5, height = 8)
  par(mfrow = c(4, 2), mar = c(4.2, 3.6, 1.6, 0.4), mgp = c(2.2, 0.6, 0), las = 2,
      cex.axis = 0.7, cex.lab = 0.85, cex.main = 0.85)
  for (s in sort(unique(out$scenario))) {
    d <- out[out$scenario == s & out$method %in% meth, ]
    d$method <- factor(d$method, levels = meth)
    boxplot(log10(PE) ~ method, data = d, names = lab, col = cols, border = "#404040",
            outline = FALSE, xlab = "", ylab = expression(log[10]~"PE"),
            main = sprintf("%s   p = %d, %s, SNR = %g", LETTERS[which(sort(unique(out$scenario)) == s)],
                           d$p[1], d$cov[1], d$snr[1]),
            boxwex = 0.6, lwd = 0.8)
  }
  dev.off()
}

###############################################################################
## Study 3: computing time
###############################################################################

study3 <- function() {
  set.seed(303)
  grid <- expand.grid(n = c(100, 200, 400), p = c(1000, 5000, 20000))
  tt <- NULL
  for (i in seq_len(nrow(grid))) {
    n <- grid$n[i]; p <- grid$p[i]
    X <- matrix(rnorm(n * p), n, p)
    beta <- numeric(p); beta[1:10] <- 1
    y <- drop(X %*% beta) + rnorm(n)
    g <- hr_grid(X, y)                       # one SVD; 12 x 40 grid
    t_hr <- system.time(hr(X, y, k = g$k, gamma = g$gamma[6, ]))[3]
    t_cvhr <- system.time(cv_hr(X, y, nfolds = 10, grid = g))[3]
    t_las <- system.time(glmnet(X, y))[3]
    t_cvlas <- system.time(cv.glmnet(X, y, nfolds = 10))[3]
    tt <- rbind(tt, data.frame(n = n, p = p, HR_path = t_hr, HR_cv = t_cvhr,
                               Lasso_path = t_las, Lasso_cv = t_cvlas))
    cat(sprintf("n=%d p=%d: HR path %.2fs, HR cv %.2fs, lasso path %.2fs, lasso cv %.2fs\n",
                n, p, t_hr, t_cvhr, t_las, t_cvlas))
  }
  write.csv(tt, "results/study3_timing.csv", row.names = FALSE)
  plot3(tt)
  tt
}

plot3 <- function(tt) {
  pdf("figures/fig3_timing.pdf", width = 7, height = 3.4)
  par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 0.5), mgp = c(2.4, 0.7, 0), las = 1,
      cex.axis = 0.85)
  for (nn in c(100, 400)) {
    d <- tt[tt$n == nn, ]
    yl <- range(d[, c("HR_cv", "Lasso_cv")])
    plot(d$p, d$HR_cv, log = "xy", type = "b", pch = 16, col = COL_HR, lwd = 2,
         ylim = yl, xlab = "p", ylab = "Seconds (10-fold CV)", bty = "l",
         main = paste(ifelse(nn == 100, "A", "B"), "  n =", nn), cex.main = 0.9)
    lines(d$p, d$Lasso_cv, type = "b", pch = 17, col = COL_GREY, lwd = 2)
    legend("topleft", c("HR", "Lasso"), col = c(COL_HR, COL_GREY), pch = c(16, 17),
           lwd = 2, bty = "n", cex = 0.8)
  }
  dev.off()
}

###############################################################################
## Study 4: block-correlated predictors
###############################################################################

## Block covariance: G blocks of equal size, correlation rho within a block
## and zero between blocks.
block_cov <- function(p, G, rho) {
  m <- p / G; S <- matrix(0, p, p)
  for (g in seq_len(G)) { idx <- ((g - 1) * m + 1):(g * m); S[idx, idx] <- rho }
  diag(S) <- 1
  S
}

## Regression vector for the block design. Block g of the active blocks has
## coefficients that sum to total[g]:
##   "constant"  all members of the block share the total equally;
##   "half"      the first half of the block shares it, the second half is zero;
##   "single"    the first member carries it, the others are zero.
block_beta <- function(p, G, active, total, pattern) {
  m <- p / G; beta <- numeric(p)
  for (a in seq_along(active)) {
    idx <- ((active[a] - 1) * m + 1):(active[a] * m)
    beta[idx] <- switch(pattern,
                        constant = total[a] / m,
                        half = c(rep(2 * total[a] / m, m / 2), rep(0, m / 2)),
                        single = c(total[a], rep(0, m - 1)))
  }
  beta
}

## Scenarios 1 to 4 are reported in the paper; scenarios 5 to 8 are the
## departures reported in the Supplementary Material: weaker within-block
## correlation, half of each active block inactive, a single active variable
## per block, and more blocks.
SCEN4 <- data.frame(p = c(200, 1000, 200, 1000, 1000, 1000, 1000, 1000),
                    snr = c(1, 1, 5, 5, 5, 5, 5, 5),
                    rho = c(0.8, 0.8, 0.8, 0.8, 0.5, 0.8, 0.8, 0.8),
                    pattern = c(rep("constant", 5), "half", "single", "constant"),
                    G = c(20, 20, 20, 20, 20, 20, 20, 50),
                    stringsAsFactors = FALSE)

study4 <- function(nrep = 100, n = 100, active = c(1, 2, 5, 10, 20),
                   total = c(3, -2, 1.5, 1, -1), scen = SCEN4) {
  res <- list()
  for (s in seq_len(nrow(scen))) {
    fn <- sprintf("results/study4_scenario%d.rds", s)
    if (file.exists(fn)) { res[[s]] <- readRDS(fn); next }
    set.seed(4000 + s)
    p <- scen$p[s]; G <- scen$G[s]
    Sigma <- block_cov(p, G, scen$rho[s])
    R <- chol(Sigma)
    beta <- block_beta(p, G, active, total, scen$pattern[s])
    sigma <- sqrt(drop(crossprod(beta, Sigma %*% beta)) / scen$snr[s])
    one <- function(r) {
      X <- rdesign(n, R); y <- drop(X %*% beta) + rnorm(n, sd = sigma)
      f <- fit_all(X, y)
      ev <- t(sapply(f$coef, evaluate, beta = beta, Sigma = Sigma))
      ## Relative distance between beta and its projection onto the row space
      ## of the centred design (no random numbers are used here).
      sv <- svd(scale(X, center = TRUE, scale = FALSE)); keep <- sv$d > 1e-10 * sv$d[1]
      theta <- drop(sv$v[, keep] %*% crossprod(sv$v[, keep], beta))
      data.frame(scenario = s, p = p, snr = scen$snr[s], rho = scen$rho[s], G = G,
                 pattern = scen$pattern[s], sigma = sigma, rep = r, method = rownames(ev), ev,
                 time = f$time[rownames(ev)],
                 reldist = sqrt(sum((theta - beta)^2) / sum(beta^2)), row.names = NULL)
    }
    t0 <- proc.time()[3]
    res[[s]] <- do.call(rbind, mclapply(seq_len(nrep), one, mc.cores = NCORES,
                                        mc.set.seed = TRUE))
    saveRDS(res[[s]], fn)
    cat(sprintf("block scenario %d/%d (p=%d, SNR=%g, rho=%g, G=%d, %s) done in %.0f s\n", s,
                nrow(scen), p, scen$snr[s], scen$rho[s], G, scen$pattern[s], proc.time()[3] - t0))
  }
  out <- do.call(rbind, res)
  saveRDS(out, "results/study4.rds")
  out
}

## Design diagnostics for Remark 4: for 50 designs per setting, the distance
## between beta and the row space of the centred design, the distance between
## beta and the span of the q = 20 leading right singular vectors, and the
## eigenvalues lambda_20, lambda_21 and lambda_r of X'X for the centred and
## scaled design used by the estimators. The sparse designs of Study 2 are
## included for contrast. No estimator is fitted here.
design4 <- function(ndes = 50, n = 100, q = 20, active = c(1, 2, 5, 10, 20),
                    total = c(3, -2, 1.5, 1, -1)) {
  set.seed(4100)
  setting <- rbind(data.frame(design = "block", p = c(200, 1000), rho = 0.8, pattern = "constant", G = 20),
                   SCEN4[5:8, c("p", "rho", "pattern", "G")] |> transform(design = "block"),
                   data.frame(design = c("AR", "AR", "CS", "CS"), p = c(200, 1000, 200, 1000),
                              rho = c(0.5, 0.5, 0.3, 0.3), pattern = "sparse", G = NA))
  tab <- do.call(rbind, lapply(seq_len(nrow(setting)), function(i) {
    p <- setting$p[i]
    if (setting$design[i] == "block") {
      Sigma <- block_cov(p, setting$G[i], setting$rho[i])
      beta <- block_beta(p, setting$G[i], active, total, setting$pattern[i])
    } else {
      Sigma <- if (setting$design[i] == "AR") ar1(p, 0.5) else cs(p, 0.3)
      beta <- numeric(p)
      beta[c(1, 2, 5, 10, 20, 50, 100, 150)] <- c(3, 1.5, 2, -2.5, 1, -1, 2, -1.5)
    }
    R <- chol(Sigma)
    d <- t(replicate(ndes, {
      X <- rdesign(n, R)
      Xc <- scale(X, center = TRUE, scale = FALSE)
      sv <- svd(Xc); keep <- sv$d > 1e-10 * sv$d[1]; V <- sv$v[, keep]
      res_row <- beta - drop(V %*% crossprod(V, beta))
      res_q <- beta - drop(V[, 1:q] %*% crossprod(V[, 1:q], beta))
      ## Eigenvalues for the centred and scaled design, as used in fitting.
      xs <- sqrt(colSums(Xc^2) / n); lam <- svd(sweep(Xc, 2, xs, "/"), nu = 0, nv = 0)$d^2
      lam <- lam[lam > 1e-20 * lam[1]]
      c(rel_row = sqrt(sum(res_row^2) / sum(beta^2)), rel_q = sqrt(sum(res_q^2) / sum(beta^2)),
        lambda_1 = lam[1], lambda_q = lam[q], lambda_q1 = lam[q + 1], lambda_r = lam[length(lam)],
        rank = length(lam))
    }))
    data.frame(setting[i, ], t(colMeans(d)))
  }))
  write.csv(tab, "results/study4_design.csv", row.names = FALSE)
  tab
}

## Mean (standard error) of every measure, and paired comparisons of the
## scaled HR estimator with each other method on the same replications.
summarise4 <- function(out, file = "results/study4_summary.csv") {
  se <- function(v) sd(v) / sqrt(length(v))
  tab <- do.call(rbind, lapply(split(out, list(out$method, out$scenario), drop = TRUE), function(d) {
    ref <- out[out$scenario == d$scenario[1] & out$method == "HR_S", ]
    dif <- d$PE[order(d$rep)] - ref$PE[order(ref$rep)]
    data.frame(scenario = d$scenario[1], p = d$p[1], snr = d$snr[1], rho = d$rho[1],
               G = d$G[1], pattern = d$pattern[1], method = d$method[1], PE = mean(d$PE), PE_se = se(d$PE),
               L2 = mean(d$L2), L2_se = se(d$L2), TPR = mean(d$TPR), FPR = mean(d$FPR),
               Size = mean(d$Size), time = mean(d$time), reldist = mean(d$reldist),
               PE_minus_HRS = mean(dif), PE_minus_HRS_se = se(dif),
               prop_HRS_better = mean(dif > 0))
  }))
  tab <- tab[order(tab$scenario, match(tab$method, METH_X)), ]
  write.csv(tab, file, row.names = FALSE)
  tab
}

## LaTeX tables for Study 4: the paper reports scenarios 1 to 4, the
## Supplementary Material scenarios 5 to 8 with the scaled hard-thresholded
## ridge added.
tables4 <- function(out) {
  hd <- function(d) sprintf("$p=%d$, %d blocks, $\\rho=%g$, %s, SNR $=%g$", d$p[1], d$G[1], d$rho[1],
                            c(constant = "constant within blocks", half = "half of each block active",
                              single = "one variable per block")[d$pattern[1]], d$snr[1])
  hd_main <- function(d) sprintf("$p=%d$, blocks of size %d, SNR $=%g$", d$p[1], d$p[1] / d$G[1], d$snr[1])
  table2(out[out$scenario %in% 1:4, ], "results/table_study4.tex", header = hd_main)
  table2(out[out$scenario %in% 5:8, ], "results/table_study4_supp.tex",
         meth = METH_X, lab = LAB_X, header = hd)
}

plot4 <- function(out, meth = METH, lab = LAB) {
  cols <- ifelse(grepl("^HR", meth), COL_HR, COL_GREY)
  pdf("figures/fig6_block_prediction_error.pdf", width = 7.5, height = 4.6)
  par(mfrow = c(2, 2), mar = c(4.2, 3.6, 1.6, 0.4), mgp = c(2.2, 0.6, 0), las = 2,
      cex.axis = 0.7, cex.lab = 0.85, cex.main = 0.85)
  for (s in 1:4) {
    d <- out[out$scenario == s & out$method %in% meth, ]
    d$method <- factor(d$method, levels = meth)
    boxplot(log10(PE) ~ method, data = d, names = lab, col = cols, border = "#404040",
            outline = FALSE, xlab = "", ylab = expression(log[10]~"PE"),
            main = sprintf("%s   p = %d, SNR = %g", LETTERS[s], d$p[1], d$snr[1]),
            boxwex = 0.6, lwd = 0.8)
  }
  dev.off()
}

## Scaled soft and scaled hard thresholding side by side in every simulation
## setting of Studies 2, 4 and 5 (Supplementary Material).
table_scaled <- function(out2, out4, out5 = NULL, file = "results/table_scaled.tex") {
  f <- function(v) sprintf("%.2f (%.2f)", mean(v), sd(v) / sqrt(length(v)))
  row <- function(d, label) {
    a <- d[d$method == "HR_S", ]; b <- d[d$method == "HardRidge_S", ]
    a <- a[order(a$rep), ]; b <- b[order(b$rep), ]
    sprintf("%s & %s & %.1f & %s & %.1f & %.2f\\\\", label, f(a$PE), mean(a$Size), f(b$PE),
            mean(b$Size), mean(a$PE < b$PE))
  }
  lines <- c()
  for (s in sort(unique(out2$scenario))) {
    d <- out2[out2$scenario == s, ]
    lines <- c(lines, row(d, sprintf("Sparse, %s, $p=%d$, SNR $=%g$",
                                     ifelse(d$cov[1] == "AR", "AR", "CS"), d$p[1], d$snr[1])))
  }
  for (s in sort(unique(out4$scenario))) {
    d <- out4[out4$scenario == s, ]
    lines <- c(lines, row(d, sprintf("Block, $G=%d$, $\\rho=%g$, %s, $p=%d$, SNR $=%g$", d$G[1],
                                     d$rho[1], d$pattern[1], d$p[1], d$snr[1])))
  }
  if (!is.null(out5)) for (s in sort(unique(out5$scenario))) {
    d <- out5[out5$scenario == s, ]
    lines <- c(lines, row(d, sprintf("%s, $p=%d$, SNR $=%g$", c(AR = "Sparse, AR", CS = "Sparse, CS",
                                     block = "Block, $G=20$, $\\rho=0.8$, constant")[d$design[1]], d$p[1], d$snr[1])))
  }
  writeLines(lines, file)
  invisible(lines)
}

###############################################################################
## Study 5: low signal-to-noise ratio
###############################################################################

## SNR = 0.2 for the four sparse designs of Study 2 (scenarios 1 to 4) and the
## two block designs of Study 4 (scenarios 5 and 6). Seeds 5000 + scenario.
SCEN5 <- data.frame(design = c("AR", "AR", "CS", "CS", "block", "block"),
                    p = c(200, 1000, 200, 1000, 200, 1000), snr = 0.2,
                    stringsAsFactors = FALSE)

study5 <- function(nrep = 100, n = 100, scen = SCEN5) {
  res <- list()
  for (s in seq_len(nrow(scen))) {
    fn <- sprintf("results/study5_scenario%d.rds", s)
    if (file.exists(fn)) { res[[s]] <- readRDS(fn); next }
    set.seed(5000 + s)
    p <- scen$p[s]
    if (scen$design[s] == "block") {
      Sigma <- block_cov(p, 20, 0.8)
      beta <- block_beta(p, 20, c(1, 2, 5, 10, 20), c(3, -2, 1.5, 1, -1), "constant")
    } else {
      Sigma <- if (scen$design[s] == "AR") ar1(p, 0.5) else cs(p, 0.3)
      beta <- numeric(p)
      beta[c(1, 2, 5, 10, 20, 50, 100, 150)] <- c(3, 1.5, 2, -2.5, 1, -1, 2, -1.5)
    }
    R <- chol(Sigma)
    signal <- drop(crossprod(beta, Sigma %*% beta))
    sigma <- sqrt(signal / scen$snr[s])
    one <- function(r) {
      X <- rdesign(n, R); y <- drop(X %*% beta) + rnorm(n, sd = sigma)
      f <- fit_all(X, y)
      ev <- t(sapply(f$coef, evaluate, beta = beta, Sigma = Sigma))
      data.frame(scenario = s, design = scen$design[s], p = p, snr = scen$snr[s],
                 signal = signal, rep = r, method = rownames(ev), ev,
                 time = f$time[rownames(ev)], row.names = NULL)
    }
    t0 <- proc.time()[3]
    res[[s]] <- do.call(rbind, mclapply(seq_len(nrep), one, mc.cores = NCORES,
                                        mc.set.seed = TRUE))
    saveRDS(res[[s]], fn)
    cat(sprintf("low-SNR scenario %d/%d (%s, p=%d) done in %.0f s\n", s, nrow(scen),
                scen$design[s], p, proc.time()[3] - t0))
  }
  out <- do.call(rbind, res)
  saveRDS(out, "results/study5.rds")
  out
}

## Supplementary table for Study 5; the null model has prediction error equal
## to the signal variance, which is given in each header.
table5 <- function(out) {
  hd <- function(d) sprintf("%s, $p=%d$, SNR $=0.2$, null-model PE $=%.2f$",
                            c(AR = "Autoregressive", CS = "Compound symmetry",
                              block = "Block")[d$design[1]], d$p[1], d$signal[1])
  table2(out, "results/table_study5.tex", meth = METH_X, lab = LAB_X, header = hd)
}

###############################################################################
## Run
###############################################################################

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  which_study <- if (length(args)) args else c("1", "2", "3", "4", "5")
  if ("1" %in% which_study) print(summarise1(study1()))
  if ("2" %in% which_study) { out <- study2(); print(summarise2(out)); plot2(out); tables2(out) }
  if ("3" %in% which_study) print(study3())
  if ("5" %in% which_study) { out5 <- study5(); table5(out5) }
  if ("4" %in% which_study) {
    out4 <- study4(); print(summarise4(out4)); plot4(out4); tables4(out4)
    print(design4())
    if (file.exists("results/study2.rds"))
      table_scaled(readRDS("results/study2.rds"), out4,
                   if (file.exists("results/study5.rds")) readRDS("results/study5.rds"))
  }
  print(sessionInfo())
}
