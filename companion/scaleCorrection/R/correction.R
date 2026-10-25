# Deterministic, balanced inner folds in a separate random-number namespace.
# The original data and outer-fold streams are never consumed or changed.
variance_folds <- function(n, unit_id, subset = "full") {
  stopifnot(length(n)==1L, n>=20L, length(unit_id)==1L, length(subset)==1L)
  key <- paste("scaleCorrection-variance-v1", unit_id, subset, sep="/")
  seed <- strtoi(substr(digest::digest(key, algo="sha256", serialize=FALSE),1,7),16L)
  existed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  if(existed) old <- get(".Random.seed", .GlobalEnv)
  kind <- RNGkind()
  on.exit({do.call(RNGkind,as.list(kind)); if(existed)assign(".Random.seed",old,.GlobalEnv)
    else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)})
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  sample(rep(seq_len(10L),length.out=n))
}

# Matches edition-4 sigma2_lasso, including the estimated intercept and floor.
# Every call uses only its supplied training observations.
estimate_variance <- function(X,y,unit_id,subset="full") {
  fid <- variance_folds(nrow(X),unit_id,subset)
  warnings <- character()
  fit <- withCallingHandlers(glmnet::cv.glmnet(X,drop(y),alpha=1,foldid=fid,
    intercept=TRUE,standardize=TRUE,parallel=FALSE),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  if(length(warnings))stop(paste("Variance CV warning:",paste(unique(warnings),collapse="; ")))
  stopifnot(fit$glmnet.fit$jerr==0L,all(is.finite(fit$cvm)))
  b <- as.numeric(stats::coef(fit,s="lambda.min")); size <- sum(b[-1]!=0)
  rss <- sum((drop(y)-b[1]-drop(X%*%b[-1]))^2)
  denominator <- max(1,nrow(X)-size-1L); sigma2 <- rss/denominator
  stopifnot(is.finite(sigma2),sigma2>=0)
  list(sigma2=sigma2,rss=rss,size=size,denominator=denominator,
    denominator_floor_used=nrow(X)-size-1L<1L,lambda=fit$lambda.min,
    foldid_sha256=digest::digest(fid,algo="sha256"),subset=subset)
}

# New joint CV surface supplies both soft-thresholded and gamma=0 controls.
# The gamma=0 control is tuned over the very same 20 ridge penalties.
corrected_pair <- function(X,y,grid,foldid,unit_id,variance_fun=estimate_variance) {
  X <- as.matrix(X); y <- drop(y); n <- nrow(X)
  softridge:::.validate_folds(foldid,n)
  stopifnot(length(grid$k)>0,all(is.finite(grid$k)),all(grid$k>0),
    length(grid$fraction)>0,grid$fraction[1]==0,all(grid$fraction>=0),all(grid$fraction<=1))
  err <- matrix(0,length(grid$k),length(grid$fraction))
  variance <- vector("list",max(foldid)+1L)
  for(f in seq_len(max(foldid))) {
    tr <- foldid!=f
    pr <- softridge:::.prep_design(X[tr,,drop=FALSE],y[tr],TRUE,TRUE)
    sv <- softridge:::.svd_design(pr$X)
    B <- softridge:::ridge_path(pr$X,pr$y,grid$k,sv)
    variance[[f]] <- variance_fun(X[tr,,drop=FALSE],y[tr],unit_id,paste0("outer-",f))
    xt <- sweep(sweep(X[!tr,,drop=FALSE],2,pr$xm,"-"),2,pr$xs,"/")
    out <- correction_surface(pr$X,pr$y,xt,y[!tr]-pr$ym,B,sv$v,sv$d,
      grid$k,grid$fraction,variance[[f]]$sigma2)
    err <- err+out$score
  }
  err <- err/n
  best <- which(err==min(err),arr.ind=TRUE)[1,]
  control <- c(which.min(err[,1]),1L)
  pr <- softridge:::.prep_design(X,y,TRUE,TRUE); sv <- softridge:::.svd_design(pr$X)
  B <- softridge:::ridge_path(pr$X,pr$y,grid$k,sv)
  variance[[length(variance)]] <- variance_fun(X,y,unit_id,"full")
  surface <- correction_surface(pr$X,pr$y,matrix(numeric(),0,ncol(X)),numeric(),
    B,sv$v,sv$d,grid$k,grid$fraction,variance[[length(variance)]]$sigma2)
  final <- function(index) {
    i<-index[1];g<-index[2];gamma<-grid$fraction[g]*max(abs(B[,i]))
    th<-softridge::soft_threshold(B[,i],gamma); cc<-surface$scale[i,g]
    beta<-cc*th/pr$xs
    list(coef=c(pr$ym-sum(beta*pr$xm),beta),k=grid$k[i],k_index=i,
      fraction=grid$fraction[g],gamma_index=g,gamma=gamma,scale=cc,
      ls_scale=surface$ls_scale[i,g],df=surface$df[i,g],cv=err[i,g],size=sum(beta!=0))
  }
  list(SSTR_c0=final(best),Ridge_c0=final(control),cv_error=err,
    variance=variance,foldid=foldid,grid=grid,projection=list(v=sv$v,xs=pr$xs),
    status="plug-in optimism correction; no adaptive unbiased-risk guarantee")
}
