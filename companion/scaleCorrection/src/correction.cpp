// Additive kernel. No changes to the frozen softridge numerical interfaces.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

// X and y are centered and RMS-standardized by the training subset only.
// The intercept contribution is excluded from df because y has been centered.
// df is the fixed-threshold expression, used here as a plug-in correction;
// response-dependent thresholds and estimated scales are not covered by SURE.
// [[Rcpp::export]]
Rcpp::List correction_surface(const arma::mat& X, const arma::vec& y,
    const arma::mat& Xtest, const arma::vec& ytest, const arma::mat& B,
    const arma::mat& V, const arma::vec& d, const arma::vec& k,
    const arma::vec& fraction, double sigma2) {
  const arma::uword nk = k.n_elem, ng = fraction.n_elem, p = X.n_cols;
  if (!nk || !ng || X.n_rows != y.n_elem || Xtest.n_rows != ytest.n_elem ||
      Xtest.n_cols != p || B.n_rows != p || B.n_cols != nk ||
      V.n_rows != p || V.n_cols != d.n_elem || !X.is_finite() ||
      !y.is_finite() || !Xtest.is_finite() || !ytest.is_finite() ||
      !B.is_finite() || !V.is_finite() || !d.is_finite() ||
      !k.is_finite() || arma::any(k <= 0) || !fraction.is_finite() ||
      arma::any(fraction < 0) || arma::any(fraction > 1) ||
      !std::isfinite(sigma2) || sigma2 < 0)
    Rcpp::stop("invalid correction inputs");
  arma::mat score(nk, ng), scale(nk, ng), df(nk, ng), ls(nk, ng);
  arma::mat V2 = arma::square(V);
  arma::vec d2 = arma::square(d);
  for (arma::uword i = 0; i < nk; ++i) {
    Rcpp::checkUserInterrupt();
    arma::vec ak = V2 * (d2 / (d2 + k[i]));
    arma::mat th(p, ng);
    const double mx = arma::abs(B.col(i)).max();
    for (arma::uword g = 0; g < ng; ++g) {
      const double threshold = fraction[g] * mx;
      double trace = 0;
      for (arma::uword j = 0; j < p; ++j) {
        const double b = B(j, i);
        th(j, g) = std::copysign(std::max(0.0, std::abs(b) - threshold), b);
        if (th(j, g) != 0) trace += ak[j];
      }
      df(i, g) = trace;
    }
    arma::mat fitted = X * th, test = Xtest * th;
    for (arma::uword g = 0; g < ng; ++g) {
      const double den = arma::dot(fitted.col(g), fitted.col(g));
      const double yfit = arma::dot(y, fitted.col(g));
      // A zero fitted vector has the same prediction for every scale.
      // Store scale zero, with no arbitrary near-zero numerical cutoff.
      const double c = den > 0 ? std::max(0.0, (yfit - sigma2 * df(i, g)) / den) : 0;
      const double cls = den > 0 ? std::max(0.0, yfit / den) : 0;
      scale(i, g) = c; ls(i, g) = cls;
      arma::vec residual = ytest - c * test.col(g);
      score(i, g) = arma::dot(residual, residual);
    }
  }
  if (!score.is_finite() || !scale.is_finite()) Rcpp::stop("nonfinite correction result");
  return Rcpp::List::create(Rcpp::Named("score") = score,
      Rcpp::Named("scale") = scale, Rcpp::Named("df") = df,
      Rcpp::Named("ls_scale") = ls);
}
