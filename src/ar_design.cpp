#include <RcppArmadillo.h>
#include <cmath>
// Internal study utilities. No estimator, tuning rule or RNG is changed.

// [[Rcpp::export(rng = false)]]
arma::mat cpp_ar_design(const arma::mat& z, double rho) {
  if (!z.n_rows || !z.n_cols || !z.is_finite() || !std::isfinite(rho) || std::abs(rho) >= 1)
    Rcpp::stop("finite innovations and abs(rho) < 1 required");
  arma::mat x = z;
  const double sd = std::sqrt(1-rho*rho);
  for (arma::uword j=1; j<z.n_cols; ++j)
    x.col(j) = rho*x.col(j-1) + sd*z.col(j);
  return x;
}

// Computes b' Sigma b for Sigma[j,k] = rho^abs(j-k), without a dense matrix.
// [[Rcpp::export(rng = false)]]
double cpp_ar_quadratic(const arma::vec& b, double rho) {
  if (!b.n_elem || !b.is_finite() || !std::isfinite(rho) || std::abs(rho) >= 1)
    Rcpp::stop("finite coefficients and abs(rho) < 1 required");
  double preceding=0, q=0;
  for (arma::uword j=0; j<b.n_elem; ++j) {
    q += b[j]*b[j] + 2*b[j]*preceding;
    preceding = rho*(preceding+b[j]);
  }
  return q;
}
