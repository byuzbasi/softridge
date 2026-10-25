test_that("native SVD respects rank and reconstructs its retained subspace", {
  X<-cbind(1:6,2*(1:6),c(1,0,1,0,1,0))
  s<-softridge:::.svd_design(X)
  expect_length(s$d,2)
  expect_equal(s$u%*%(s$d*t(s$v)),X,tolerance=1e-10)
})
test_that("native CV refits keep the column-order tie rule", {
  X<-cbind(c(-2,-1,0,1,2),c(1,-1,0,-1,1))
  y<-X[,1];B<-matrix(c(1,1),2,1)
  e<-softridge:::cpp_cv_scores(X,y,X,y,B,matrix(0,1,0),1L)
  expect_equal(as.numeric(e),0,tolerance=1e-12)
})
test_that("native thresholding retains names, dimensions and missing-value types", {
  x<-matrix(c(-2,0,1,NA_real_),2,dimnames=list(c('a','b'),c('x','y')))
  expect_identical(dimnames(soft_threshold(x,.5)),dimnames(x))
  expect_equal(soft_threshold(x,.5),sign(x)*pmax(abs(x)-.5,0))
  expect_true(is.nan(soft_threshold(NaN,.5)))
  expect_identical(soft_threshold(numeric(),.5),numeric())
  expect_error(softridge:::cpp_soft(1,numeric()),"nonnegative")
  X<-matrix(1:12,4);y<-1:4;g<-hr_grid(X,y,nk=2,ng=2);g$gamma[]<--1
  expect_error(sure_hr(X,y,1,g),"nonnegative")
})
