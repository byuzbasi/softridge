test_that("rank-deficient refitting returns the minimum-norm coefficients", {
  X<-cbind(c(-2,-1,0,1,2),c(-2,-1,0,1,2));y<-X[,1]
  f<-softridge:::.mp_refit(X,y)
  expect_equal(f$coef,c(0.5,0.5),tolerance=1e-10)
  expect_equal(f$rank,1L)
  expect_equal(drop(X%*%f$coef),y,tolerance=1e-10)
})
test_that("fraction grids cannot import full-response absolute thresholds into folds", {
  set.seed(41);X<-matrix(rnorm(24*30),24,30);y<-2*X[,1]+rnorm(24)
  folds<-rep(1:4,6);g<-hr_grid(X,y,nk=2,ng=4)
  a<-cv_hr(X,y,grid=g,foldid=folds)
  g$gamma[]<-1e8
  b<-cv_hr(X,y,grid=g,foldid=folds)
  expect_equal(a$cv_error,b$cv_error)
  expect_equal(a$coef,b$coef)
})
test_that("refit candidate size respects the smallest training rank", {
  set.seed(42);X<-matrix(rnorm(16*20),16,20);y<-rnorm(16)
  f<-cv_hr(X,y,grid=hr_grid(X,y,nk=2,ng=3),foldid=rep(1:4,4),refit=TRUE,max_size=40)
  expect_lte(f$refit_diagnostics$effective_max_size,11)
  expect_true(all(is.finite(f$coef)))
  expect_error(cv_hr(X,y,foldid=c(rep(1,8),rep(3,8))),"consecutive")
})
test_that("fixed-tuning SURE divergence equals a numerical Jacobian trace", {
  set.seed(43);X<-matrix(rnorm(12*18),12,18);y<-rnorm(12);h<-1e-6
  fit<-function(z)predict(hr(X,z,k=2,gamma=.07),X)
  trace<-sum(vapply(seq_along(y),function(i) {
    yp<-ym<-y;yp[i]<-yp[i]+h;ym[i]<-ym[i]-h
    (fit(yp)[i]-fit(ym)[i])/(2*h)
  },numeric(1)))
  expect_equal(as.numeric(hr_sure(X,y,1,2,.07)$df),trace,tolerance=1e-6)
})
test_that("right risk derivative agrees with the stated local criterion", {
  set.seed(44);X<-matrix(rnorm(9*13),9,13);b<-c(1,-1,rep(0,11));k<-2
  q<-hr_risk(X,b,1,k,0)
  delta<-sum((q$m-q$theta)*(2*pnorm(q$m/q$s)-1)+2*q$s*dnorm(q$m/q$s))
  numerical<-(hr_risk(X,b,1,k,1e-6)$risk-q$risk)/1e-6
  expect_equal(numerical,-2*delta,tolerance=1e-4)
})
test_that("single predictor shapes and invalid thresholds are handled", {
  X<-matrix(seq_len(10),ncol=1);y<-seq_len(10)+sin(seq_len(10))
  f<-hr(X,y,k=1,gamma=.1)
  expect_length(coef(f),2)
  g<-hr_grid(X,y,nk=2,ng=1)
  expect_equal(dim(g$gamma),c(2L,1L))
  expect_length(cv_hr(X,y,grid=g,foldid=rep(1:2,5))$coef,2L)
  expect_length(sure_hr(X,y,1,grid=g)$coef,2L)
  expect_error(hr(X,y,k=1,gamma=-1,nonneg=TRUE),"nonnegative")
})

test_that("theoretical routines reject parameters outside their assumptions", {
  X<-diag(2);b<-c(1,1)
  expect_error(hr_risk(X,b,1,-1,.1),"positive")
  expect_error(hr_risk(X,b,0,1,.1),"positive")
  expect_error(hr_risk(X,b,1,1,-.1),"nonnegative")
  expect_error(hr_risk(X,1,1,1,.1),"matching")
  expect_error(hr_sure(X,b,-1,1,.1),"positive")
  expect_error(sure_hr(X,b,NA_real_),"positive")
  expect_error(hr_grid(X,b,nk=0),"invalid grid")
})
