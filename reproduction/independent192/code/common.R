# Comparator fitting functions for the independent study.
# Estimator formulas, CV grids and convergence rules are retained.
suppressPackageStartupMessages({library(softridge); library(glmnet); library(ncvreg)})

## Hard-thresholded ridge estimator of Shao and Deng (2012), tuned by
## K-fold cross-validation over the same (k, gamma) grid as HR, so that the
## only difference between the two is the thresholding rule.
cv_hard_ridge <- function(X,y,grid,foldid,scale=FALSE) {
  f <- softridge:::.cv_hard_ridge(X,y,grid,foldid,scale)
  cf <- f$coef; attr(cf,"tuning") <- f; cf
}

fit_all <- function(X, y, nfolds = 10, foldid = NULL, grid = NULL, max_size = 40) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(foldid)) foldid <- sample(rep(seq_len(nfolds), length.out = n))
  softridge:::.validate_folds(foldid, n)
  out <- list(); tm <- c(); diagnostics <- list()
  tic <- function() proc.time()[3]

  t0 <- tic(); f <- cv.glmnet(X, y, alpha = 0, foldid = foldid, lambda.min.ratio=1e-6, nlambda=200)
  out$Ridge <- as.numeric(coef(f, s = "lambda.min")); tm["Ridge"] <- tic() - t0

  t0 <- tic(); flas <- cv.glmnet(X, y, alpha = 1, foldid = foldid)
  out$Lasso <- as.numeric(coef(flas, s = "lambda.min")); tm["Lasso"] <- tic() - t0

  ## Elastic net: the mixing parameter is chosen on the same folds.
  t0 <- tic(); cvs <- lapply(c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8), function(a)
    cv.glmnet(X, y, alpha = a, foldid = foldid))
  j <- which.min(sapply(cvs, function(f) min(f$cvm)))
  out$ENet <- as.numeric(coef(cvs[[j]], s = "lambda.min")); tm["ENet"] <- tic() - t0

  t0 <- tic(); fs <- cv_ncvreg_audited(X, y, penalty = "SCAD", foldid = foldid)
  out$SCAD <- fs$coef; diagnostics$SCAD <- fs$diagnostics; tm["SCAD"] <- tic() - t0
  t0 <- tic(); fm <- cv_ncvreg_audited(X, y, penalty = "MCP", foldid = foldid)
  out$MCP <- fm$coef; diagnostics$MCP <- fm$diagnostics; tm["MCP"] <- tic() - t0

  if (is.null(grid)) grid <- hr_grid(X, y)
  t0 <- tic(); out$HardRidge <- cv_hard_ridge(X, y, grid, foldid); tm["HardRidge"] <- tic() - t0
  t0 <- tic(); cvh <- cv_hr(X, y, grid = grid, foldid = foldid)
  diagnostics$HR <- cvh
  out$HR <- cvh$coef; tm["HR"] <- tic() - t0
  diagnostics$HardRidge <- attr(out$HardRidge,"tuning")
  t0 <- tic(); sc <- cv_hr(X,y,grid=grid,foldid=foldid,scale=TRUE)
  out$HR_S <- sc$coef; diagnostics$HR_S <- sc; tm["HR_S"] <- tic()-t0
  t0 <- tic(); hs <- cv_hard_ridge(X,y,grid,foldid,scale=TRUE)
  out$HardRidge_S <- hs; diagnostics$HardRidge_S <- attr(hs,"tuning"); tm["HardRidge_S"] <- tic()-t0
  diagnostics$glmnet_surfaces <- lapply(list(Ridge=f,Lasso=flas,ENet=cvs[[j]]),function(z)
    list(lambda=z$lambda,cvm=z$cvm,cvsd=z$cvsd,lambda_min=z$lambda.min,nzero=z$nzero))
  diagnostics$enet <- list(alpha_grid=c(0.01,0.05,0.1,0.2,0.5,0.8),selected_alpha=c(0.01,0.05,0.1,0.2,0.5,0.8)[j],
    alpha_minimum=vapply(cvs,function(z)min(z$cvm),numeric(1)))
  diagnostics$glmnet <- list(ridge_jerr = f$glmnet.fit$jerr, lasso_jerr = flas$glmnet.fit$jerr,
    enet_jerr = vapply(cvs, function(z) z$glmnet.fit$jerr, numeric(1)))
  status <- setNames(rep("ok", length(out)), names(out))
  order <- c("Ridge","Lasso","ENet","SCAD","MCP","HardRidge","HR","HR_S","HardRidge_S")
  list(coef = out[order], time = tm[order], foldid = foldid, diagnostics = diagnostics, status = status[order])
}


capture_solver<-function(expr) {
 w<-character();value<-withCallingHandlers(expr,warning=function(x){w<<-c(w,conditionMessage(x));invokeRestart('muffleWarning')})
 if(any(!grepl('Maximum number of iterations reached',w,fixed=TRUE)))stop(paste('Unexpected warning:',paste(w,collapse='; ')))
 list(value=value,warnings=w)
}
# Calls in unchanged fit_all resolve to this audited version in this script only.
cv_ncvreg_audited<-function(X,y,penalty,foldid) {
 one<-function(...)capture_solver(ncvreg(...,max.iter=100000L,eps=1e-4))
 full<-one(X,y,penalty=penalty);fit<-full$value;lambda<-fit$lambda;L<-length(lambda)
 E<-matrix(NA_real_,length(y),L);iterations<-matrix(NA_integer_,max(foldid)+1L,L);iterations[1,]<-fit$iter
 warnings<-full$warnings
 for(f in seq_len(max(foldid))) {
  tr<-foldid!=f;ff<-one(X[tr,,drop=FALSE],y[tr],penalty=penalty,lambda=lambda,returnX=FALSE,convex=FALSE)
  a<-ff$value;nl<-length(a$lambda);stopifnot(max(abs(a$lambda-lambda[seq_len(nl)]))<1e-12)
  E[!tr,seq_len(nl)]<-(y[!tr]-matrix(predict(a,X[!tr,,drop=FALSE]),nrow=sum(!tr)))^2
  iterations[f+1,seq_len(nl)]<-a$iter;warnings<-c(warnings,ff$warnings)
 }
 available<-which(colSums(is.finite(E))==length(y));err<-rep(NA_real_,L)
 err[available]<-colMeans(E[,available,drop=FALSE]);j<-if(length(available))available[which.min(err[available])] else NA_integer_
 complete<-length(available)==L&&all(rowSums(iterations,na.rm=TRUE)<100000)&&!length(warnings)
 list(coef=if(is.na(j))rep(NA_real_,ncol(X)+1) else as.numeric(coef(fit,which=j)),
 diagnostics=list(selected_lambda=if(is.na(j))NA_real_ else lambda[j],index=j,lambda=lambda,iterations=iterations,
 omitted_path_indices=setdiff(seq_len(L),available),cv_error=err,squared_errors=E,
 selected_converged=complete,warnings=warnings,full_coef=as.matrix(coef(fit))))
}

source("code/oracles.R")
