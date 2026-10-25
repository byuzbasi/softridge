# Supplementary simulation only. Truth never enters feasible model fitting.
support_metrics <- function(b, beta) {
  if (anyNA(b)) return(c(TP=NA_real_,TN=NA_real_,FP=NA_real_,FN=NA_real_,MCC=NA_real_))
  selected <- b != 0; active <- beta != 0
  tp <- sum(selected & active); tn <- sum(!selected & !active)
  fp <- sum(selected & !active); fn <- sum(!selected & active)
  den <- sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  c(TP=tp,TN=tn,FP=fp,FN=fn,MCC=if(den==0)NA_real_ else mltools::mcc(TP=tp,FP=fp,TN=tn,FN=fn))
}
add_oracles <- function(fit, X, y, beta, Sigma, grid, policy) {
  stopifnot(identical(policy,"original_full_rank_OLS_otherwise_not_applicable"))
  p <- ncol(X); active <- which(beta!=0)
  start <- proc.time()[3]
  design <- cbind(1,X[,active,drop=FALSE]); q <- qr(design)
  available <- q$rank==length(active)+1L
  cf <- rep(NA_real_,p+1L)
  if(available){cfs <- qr.coef(q,y); cf <- numeric(p+1L); cf[1] <- cfs[1]; cf[active+1L] <- cfs[-1]}
  fit$coef$Oracle_OLS <- cf
  fit$time["Oracle_OLS"] <- unname(proc.time()[3]-start)
  fit$diagnostics$Oracle_OLS <- list(active=active,rank=q$rank,required_rank=length(active)+1L,available=available,
    reason=if(available) "full_rank_original_OLS" else "support_plus_intercept_design_not_full_column_rank")
  fit$status <- c(fit$status,Oracle_OLS=if(available) "ok" else "not_applicable_singular_support")
  fit
}
