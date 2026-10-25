###############################################################################
## realdata.R
##
## Real-data analysis for "High-dimensional ridge shrinkage and selection by
## soft thresholding" (Arashi and Yuzbasi): riboflavin production by
## Bacillus subtilis.
##
## Data source: the riboflavin data of Buhlmann, Kalisch and Meier (2014),
##   distributed with the R package hdi (object `riboflavin`, file
##   data/riboflavin.RData in the hdi source). n = 71 samples, p = 4088
##   log gene-expression levels, response = log riboflavin production rate.
##   Set RIBO_PATH below to the location of riboflavin.RData, or install hdi
##   and use data(riboflavin, package = "hdi").
##
## The 71 samples come from 28 fermentation batches (identifier "Fbat..." in
## the sample names), each measured at two to four time points. Samples from
## the same batch are not independent, so every split below, both the
## train/test splits and the cross-validation folds, keeps whole batches
## together.
##
## Requirements: R (>= 4.1), glmnet, ncvreg, softridge (>= 0.3.0).
## Output: results/ribo_*.csv|rds and figures/fig4_*.pdf, fig5_*.pdf.
###############################################################################

suppressPackageStartupMessages({
  library(softridge); library(glmnet); library(ncvreg); library(parallel)
})
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
NCORES <- max(1, min(2, detectCores()))
RNGkind("L'Ecuyer-CMRG")
set.seed(404)
COL_HR <- "#2a78d6"; COL_GREY <- "#8c8c8c"

RIBO_PATH <- "../data/riboflavin.RData"
if (file.exists(RIBO_PATH)) load(RIBO_PATH) else data(riboflavin, package = "hdi")
X <- as.matrix(unclass(riboflavin$x)); y <- riboflavin$y
batch <- sub(".*(Fbat[0-9]+).*", "\\1", rownames(X))
stopifnot(nrow(X) == 71, ncol(X) == 4088, length(unique(batch)) == 28)

###############################################################################
## Data profile (reported in the paper)
###############################################################################
profile <- list(n = nrow(X), p = ncol(X), missing = sum(is.na(X)) + sum(is.na(y)),
                n_batches = length(unique(batch)),
                y_summary = summary(y), y_sd = sd(y),
                shapiro_y = shapiro.test(y)$p.value,
                rank_centred = qr(scale(X, scale = FALSE))$rank,
                min_col_sd = min(apply(X, 2, sd)))
## Eigenvalues of X'X for the centred and scaled design: how much of their
## sum the leading directions carry (Remark 4 of the paper).
eig <- local({
  Xc <- scale(X, center = TRUE, scale = FALSE); xs <- sqrt(colSums(Xc^2) / nrow(X))
  lam <- svd(sweep(Xc, 2, xs, "/"), nu = 0, nv = 0)$d^2
  lam[lam > 1e-10 * lam[1]]
})
profile$eig_leading <- eig[1:5]
profile$eig_smallest <- eig[length(eig)]
profile$eig_share <- cumsum(eig)[c(1, 3, 5, 10)] / sum(eig)
print(profile)
saveRDS(profile, "results/ribo_profile.rds")

###############################################################################
## Methods (identical to those of simulation.R, cross-validation folds given)
###############################################################################
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
fit_all <- function(X, y, foldid) {
  n <- nrow(X); p <- ncol(X)
  out <- list(); tm <- c()
  tic <- function() proc.time()[3]

  ## Ridge regression over a long penalty path (six orders of magnitude), as
  ## in simulation.R.
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


## Batch-wise fold assignment: whole batches are allocated to folds.
batch_folds <- function(b, nfolds) {
  ub <- unique(b); fb <- sample(rep(seq_len(nfolds), length.out = length(ub)))
  fb[match(b, ub)]
}

###############################################################################
## Repeated batch-wise train/test splits
###############################################################################
NSPLIT <- 100
ub <- unique(batch)
one_split <- function(r) {
  test_b <- sample(ub, 7)                      # 7 of 28 batches (25%) held out
  te <- batch %in% test_b; tr <- !te
  fid <- batch_folds(batch[tr], 7)             # 7-fold batch-wise CV
  f <- fit_all(X[tr, ], y[tr], fid)
  mse <- sapply(f$coef, function(cf) mean((y[te] - cf[1] - X[te, ] %*% cf[-1])^2))
  size <- sapply(f$coef, function(cf) sum(cf[-1] != 0))
  sel <- lapply(f$coef, function(cf) which(cf[-1] != 0))
  list(res = data.frame(split = r, method = names(mse), mse = mse, size = size,
                        time = f$time[names(mse)], row.names = NULL), sel = sel)
}
t0 <- proc.time()[3]
splits <- mclapply(seq_len(NSPLIT), one_split, mc.cores = NCORES, mc.set.seed = TRUE)
cat(sprintf("splits done in %.0f s\n", proc.time()[3] - t0))
res <- do.call(rbind, lapply(splits, `[[`, "res"))
saveRDS(splits, "results/ribo_splits.rds")
tab <- aggregate(cbind(mse, size, time) ~ method, data = res, FUN = mean)
tab_se <- aggregate(cbind(mse, size) ~ method, data = res,
                    FUN = function(v) sd(v) / sqrt(length(v)))
names(tab_se)[-1] <- paste0(names(tab_se)[-1], "_se")
tab <- merge(tab, tab_se)
tab$mse_median <- aggregate(mse ~ method, data = res, FUN = median)$mse
tab <- tab[order(tab$mse), ]
write.csv(tab, "results/ribo_summary.csv", row.names = FALSE)
print(tab)

## Paired comparison of HR-OLS with each competitor on the same splits:
## proportion of splits where HR-OLS has the smaller test error.
w <- reshape(res[, c("split", "method", "mse")], idvar = "split",
             timevar = "method", direction = "wide")
names(w) <- sub("mse.", "", names(w))
paired <- sapply(setdiff(unique(res$method), "HR_OLS"),
                 function(m) mean(w$HR_OLS < w[[m]]))
print(round(paired, 2))
write.csv(data.frame(method = names(paired), prop_HR_OLS_better = paired),
          "results/ribo_paired.csv", row.names = FALSE)

## Paired difference between HR and ridge regression (reported in the paper).
d_hr <- w$HR - w$Ridge
cat(sprintf("HR minus Ridge: %.4f (%.4f)\n", mean(d_hr), sd(d_hr) / sqrt(length(d_hr))))
cat(sprintf("HR-OLS better than the lasso in %d splits\n", sum(w$HR_OLS < w$Lasso)))

## Paired comparison of the scaled HR estimator with each other method: mean
## difference in test error (method minus HR-S), its standard error and the
## proportion of splits where HR-S has the smaller test error.
paired_s <- t(sapply(setdiff(unique(res$method), "HR_S"), function(m) {
  d <- w[[m]] - w$HR_S
  c(diff = mean(d), se = sd(d) / sqrt(length(d)), prop_HR_S_better = mean(d > 0))
}))
print(round(paired_s, 4))
write.csv(data.frame(method = rownames(paired_s), paired_s), "results/ribo_paired_scaled.csv",
          row.names = FALSE)

## Selection frequencies over the splits.
selfreq <- function(m) tabulate(unlist(lapply(splits, function(s) s$sel[[m]])),
                                nbins = ncol(X)) / NSPLIT
freq <- sapply(c("HR_OLS", "Lasso", "MCP", "SCAD", "BAR"), selfreq)
rownames(freq) <- colnames(X)
top <- order(freq[, "HR_OLS"], decreasing = TRUE)[1:15]
write.csv(data.frame(gene = colnames(X)[top], freq[top, ]), "results/ribo_selfreq.csv",
          row.names = FALSE)

###############################################################################
## Full-data fits
###############################################################################
fid_full <- batch_folds(batch, 7)
full <- fit_all(X, y, fid_full)
full_sel <- lapply(full$coef, function(cf) colnames(X)[cf[-1] != 0])
saveRDS(full, "results/ribo_full.rds")
cat("\nGenes selected on the full data:\n")
for (m in c("HR_OLS", "HR_S", "Lasso", "SCAD", "MCP", "BAR", "ALasso"))
  cat(sprintf("%-7s (%2d): %s\n", m, length(full_sel[[m]]),
              paste(head(full_sel[[m]], 20), collapse = " ")))
ov <- sapply(full_sel, function(s) length(intersect(s, full_sel$HR_OLS)))
print(ov)
cvfull <- cv_hr(X, y, foldid = fid_full, refit = TRUE, max_size = 40)
cat(sprintf("HR-OLS full data: k = %.3g, model size = %d\n", cvfull$k, cvfull$size))
fl_full <- cv.glmnet(X, y, alpha = 1, foldid = fid_full)
cvscaled <- cv_hr(X, y, foldid = fid_full, scale = TRUE,
                  sigma2 = fold_variance(X, y, fid_full, sigma2_cvlasso))
cat(sprintf("HR-S full data: k = %.3g, gamma = %.3g, scale = %.3f, model size = %d\n",
            cvscaled$k, cvscaled$gamma, cvscaled$scale, cvscaled$size))

###############################################################################
## Figures
###############################################################################
meth <- c("Ridge", "Lasso", "ENet", "ALasso", "SCAD", "MCP", "BAR",
          "HardRidge", "HR", "HR_SURE", "HR_OLS", "Ridge_S", "HR_S")
lab <- c("Ridge", "Lasso", "ENet", "ALasso", "SCAD", "MCP", "BAR",
         "TRidge", "HR", "HR-SURE", "HR-OLS", "Ridge-S", "HR-S")
pdf("figures/fig4_ribo_mse.pdf", width = 7, height = 3.6)
par(mar = c(5, 4, 0.8, 0.5), mgp = c(2.5, 0.7, 0), las = 2, cex.axis = 0.8)
d <- res[res$method %in% meth, ]; d$method <- factor(d$method, levels = meth)
boxplot(mse ~ method, data = d, names = lab, xlab = "", ylab = "Test mean squared error",
        col = ifelse(grepl("^HR", meth), COL_HR, COL_GREY), border = "#404040",
        outline = FALSE, boxwex = 0.6, lwd = 0.8)
dev.off()

pdf("figures/fig5_ribo_selection.pdf", width = 7, height = 3.8)
par(mar = c(6.2, 4, 0.8, 0.5), mgp = c(2.5, 0.7, 0), las = 2, cex.axis = 0.75)
barplot(t(freq[top, c("HR_OLS", "Lasso", "MCP")]), beside = TRUE,
        col = c(COL_HR, "#eb6834", COL_GREY), border = NA, ylim = c(0, 1),
        names.arg = colnames(X)[top], ylab = "Selection frequency")
legend("topright", c("HR-OLS", "Lasso", "MCP"), fill = c(COL_HR, "#eb6834", COL_GREY),
       border = NA, bty = "n", cex = 0.8)
dev.off()

print(sessionInfo())
## LaTeX table for the paper: test error, model size and paired comparison.
tab <- read.csv("results/ribo_summary.csv"); pr <- read.csv("results/ribo_paired.csv")
lab <- c(Ridge = "Ridge", Lasso = "Lasso", ENet = "ENet", ALasso = "ALasso", SCAD = "SCAD",
         MCP = "MCP", BAR = "BAR", HardRidge = "TRidge", HR = "HR", HR_SURE = "HR-SURE",
         HR_OLS = "HR-OLS", HR_S = "HR-S", HardRidge_S = "TRidge-S",
         Ridge_S = "Ridge-S", HR_LS = "HR-LS")
tab$better <- pr$prop_HR_OLS_better[match(tab$method, pr$method)]
lines <- sprintf("%s & %.3f (%.3f) & %.3f & %.1f & %s\\\\", lab[tab$method], tab$mse, tab$mse_se,
                 tab$mse_median, tab$size,
                 ifelse(is.na(tab$better), "--", sprintf("%.2f", tab$better)))
writeLines(lines, "results/table_ribo.tex")
