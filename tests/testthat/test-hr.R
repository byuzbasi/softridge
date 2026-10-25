test_that("ridge path matches the p-dimensional formula", {
  set.seed(1)
  X <- matrix(rnorm(20 * 40), 20, 40); y <- rnorm(20); k <- 3
  b1 <- drop(softridge:::ridge_path(X, y, k))
  b2 <- drop(solve(crossprod(X) + k * diag(40), crossprod(X, y)))
  expect_equal(b1, b2, tolerance = 1e-8)
})

test_that("soft thresholding is correct", {
  expect_equal(soft_threshold(c(-2, -0.5, 0.3, 1.7), 1), c(-1, 0, 0, 0.7))
})

test_that("exact risk agrees with independent one-dimensional quadrature", {
  X<-rbind(c(1,2,0),c(0,1,2));beta<-c(1,-1,0.5);k<-2;gamma<-.2
  H<-solve(crossprod(X)+k*diag(3),t(X));mu<-drop(H%*%X%*%beta)
  sdj<-sqrt(rowSums(H^2));theta<-drop(t(X)%*%solve(tcrossprod(X),X%*%beta))
  numerical<-sum(vapply(1:3,function(j)integrate(function(z)
    (sign(z)*pmax(abs(z)-gamma,0)-theta[j])^2*dnorm(z,mu[j],sdj[j]),
    -Inf,Inf,rel.tol=1e-8,subdivisions=1000)$value,numeric(1)))
  expect_equal(hr_risk(X,beta,1,k,gamma)$risk,numerical,tolerance=1e-6)
})

test_that("hr returns coefficients that reproduce predictions", {
  set.seed(3)
  X <- matrix(rnorm(30 * 60), 30, 60); y <- X[, 1] * 2 + rnorm(30)
  f <- hr(X, y, k = 5, gamma = 0.05)
  expect_equal(length(coef(f)), 61)
  expect_equal(predict(f, X), drop(cbind(1, X) %*% coef(f)))
})

test_that("cv_hr and sure_hr run", {
  set.seed(4)
  X <- matrix(rnorm(40 * 80), 40, 80); y <- X[, 1:2] %*% c(2, -2) + rnorm(40)
  cv <- cv_hr(X, y, nfolds = 4)
  expect_true(length(cv$coef) == 81)
  cvr <- cv_hr(X, y, nfolds = 4, refit = TRUE, max_size = 10)
  expect_true(cvr$size <= 10)
  s <- sure_hr(X, y, sigma2 = 1)
  expect_true(is.finite(s$sure))
})
