#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {
Rcpp::NumericVector vec_r(const arma::vec& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}
struct ThinSVD {
  arma::mat u, v;
  arma::vec d;
};
ThinSVD decompose(const arma::mat& x, double tol = 1e-10) {
  if (x.n_rows == 0 || x.n_cols == 0 || !x.is_finite())
    Rcpp::stop("finite nonempty design required");
  ThinSVD s;
  if (!arma::svd_econ(s.u, s.d, s.v, x, "both", "dc"))
    Rcpp::stop("singular value decomposition failed");
  arma::uvec keep = arma::find(s.d > tol * s.d.max());
  s.u = s.u.cols(keep); s.v = s.v.cols(keep); s.d = s.d.elem(keep);
  return s;
}
arma::vec threshold(const arma::vec& b, double g, bool nonneg = false) {
  arma::vec out(b.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < b.n_elem; ++j) {
    if (nonneg) out[j] = std::max(b[j] - g, 0.0);
    else if (b[j] > g) out[j] = b[j] - g;
    else if (b[j] < -g) out[j] = b[j] + g;
  }
  return out;
}
double ls_scale(const arma::mat& x, const arma::vec& y, const arma::vec& th) {
  arma::vec fit = x * th;
  double den = arma::dot(fit, fit);
  if (den == 0.0) return 1.0;
  double c = std::max(0.0, arma::dot(fit, y) / den);
  if (!std::isfinite(c) || !std::isfinite(den)) Rcpp::stop("non-finite least-squares scale");
  return c;
}
arma::vec threshold_rule(const arma::vec& b, double g, bool hard) {
  if (!hard) return threshold(b, g);
  arma::vec th = b;
  th.elem(arma::find(arma::abs(b) <= g)).zeros();
  return th;
}
arma::vec refit(const arma::mat& x, const arma::vec& y, unsigned int& rank) {
  ThinSVD s = decompose(x);
  rank = s.d.n_elem;
  if (!rank) Rcpp::stop("selected design has zero numerical rank");
  arma::vec out = s.v * ((s.u.t() * y) / s.d);
  if (!out.is_finite()) Rcpp::stop("non-finite least-squares refit");
  return out;
}
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_svd(const arma::mat& x, double tol) {
  ThinSVD s = decompose(x, tol);
  return Rcpp::List::create(Rcpp::Named("u")=s.u, Rcpp::Named("d")=vec_r(s.d),
                            Rcpp::Named("v")=s.v);
}

// [[Rcpp::export(rng = false)]]
arma::mat cpp_ridge(const arma::mat& u, const arma::vec& d, const arma::mat& v,
                    const arma::vec& y, const arma::vec& k) {
  if (u.n_rows != y.n_elem || u.n_cols != d.n_elem || v.n_cols != d.n_elem)
    Rcpp::stop("incompatible SVD dimensions");
  if (!k.n_elem || !k.is_finite() || arma::any(k <= 0)) Rcpp::stop("k must be positive and finite");
  arma::vec uty = u.t() * y;
  arma::mat out(v.n_rows, k.n_elem);
  for (arma::uword i = 0; i < k.n_elem; ++i) {
    Rcpp::checkUserInterrupt();
    out.col(i) = v * ((d / (arma::square(d) + k[i])) % uty);
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector cpp_soft(Rcpp::NumericVector z, Rcpp::NumericVector gamma) {
  const R_xlen_t nz = z.size(), ng = gamma.size();
  if (!ng) Rcpp::stop("gamma must be finite and nonnegative");
  for (R_xlen_t i = 0; i < ng; ++i)
    if (!R_FINITE(gamma[i]) || gamma[i] < 0) Rcpp::stop("gamma must be finite and nonnegative");
  if (!nz) return Rcpp::clone(z);
  const R_xlen_t n = std::max(nz, ng);
  if (n % std::min(nz, ng)) Rcpp::warning("longer object length is not a multiple of shorter object length");
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    double b = z[i % nz], g = gamma[i % ng];
    if (ISNAN(b)) out[i] = b;
    else if (b > g) out[i] = b - g;
    else if (b < -g) out[i] = b + g;
    else out[i] = 0;
  }
  if (n == nz) SHALLOW_DUPLICATE_ATTRIB(out, z);
  else if (n == ng) SHALLOW_DUPLICATE_ATTRIB(out, gamma);
  return out;
}

// [[Rcpp::export(rng = false)]]
arma::vec cpp_ls_scale(const arma::mat& x, const arma::vec& y, const arma::mat& th) {
  if (x.n_rows != y.n_elem || x.n_cols != th.n_rows) Rcpp::stop("incompatible scale dimensions");
  arma::vec out(th.n_cols);
  for (arma::uword j = 0; j < th.n_cols; ++j) out[j] = ls_scale(x, y, th.col(j));
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_hr_path(const arma::mat& b, const arma::vec& gamma,
                       const arma::vec& xs, const arma::vec& xm, double ym, bool nonneg,
                       const arma::mat& x, const arma::vec& y, bool scale) {
  if (xs.n_elem != b.n_rows || xm.n_elem != b.n_rows || arma::any(xs == 0))
    Rcpp::stop("incompatible scaling dimensions");
  Rcpp::NumericVector beta(b.n_rows * b.n_cols * gamma.n_elem);
  beta.attr("dim") = Rcpp::IntegerVector::create(b.n_rows, b.n_cols, gamma.n_elem);
  arma::mat a0(b.n_cols, gamma.n_elem), sc(b.n_cols, gamma.n_elem, arma::fill::ones);
  for (arma::uword g = 0; g < gamma.n_elem; ++g) {
    Rcpp::checkUserInterrupt();
    for (arma::uword i = 0; i < b.n_cols; ++i) {
      arma::vec th = threshold(b.col(i), gamma[g], nonneg);
      if (scale) { sc(i,g) = ls_scale(x,y,th); th *= sc(i,g); }
      th /= xs;
      std::copy(th.begin(), th.end(), beta.begin() + b.n_rows * (i + b.n_cols * g));
      a0(i, g) = ym - arma::dot(th, xm);
    }
  }
  return Rcpp::List::create(Rcpp::Named("beta")=beta, Rcpp::Named("a0")=a0, Rcpp::Named("scale")=sc);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_risk(const arma::mat& x, const arma::vec& beta, double sigma,
                    double k, double gamma, bool nonneg, double scale = 1.0) {
  if (x.n_cols != beta.n_elem) Rcpp::stop("finite design and matching beta required");
  ThinSVD sv = decompose(x);
  arma::vec d2 = arma::square(sv.d);
  arma::vec theta = sv.v * (sv.v.t() * beta);
  arma::vec m = sv.v * ((d2 / (d2 + k)) % (sv.v.t() * theta));
  arma::vec s = sigma * arma::sqrt(arma::square(sv.v) * (d2 / arma::square(d2 + k)));
  if (arma::any(s <= 0)) Rcpp::stop("zero column in X: coordinatewise variance vanishes");
  arma::vec mean(m.n_elem), second(m.n_elem), mse(m.n_elem);
  for (arma::uword j = 0; j < m.n_elem; ++j) {
    double mu1 = m[j] - gamma, mu2 = -m[j] - gamma;
    double p1 = R::pnorm(mu1 / s[j], 0, 1, true, false);
    double d1 = R::dnorm(mu1 / s[j], 0, 1, false);
    mean[j] = mu1 * p1 + s[j] * d1;
    second[j] = (mu1 * mu1 + s[j] * s[j]) * p1 + mu1 * s[j] * d1;
    if (!nonneg) {
      double p2 = R::pnorm(mu2 / s[j], 0, 1, true, false);
      double d2j = R::dnorm(mu2 / s[j], 0, 1, false);
      mean[j] -= mu2 * p2 + s[j] * d2j;
      second[j] += (mu2 * mu2 + s[j] * s[j]) * p2 + mu2 * s[j] * d2j;
    }
    mse[j] = second[j] - 2 * theta[j] * mean[j] + theta[j] * theta[j];
  }
  double kappa1 = arma::dot(theta,mean), kappa2 = arma::accu(second);
  double copt = kappa2 > 0 ? std::max(kappa1/kappa2,0.0) : NA_REAL;
  double ropt = kappa2 > 0 ? arma::dot(theta,theta)-std::pow(std::max(kappa1,0.0),2)/kappa2 : NA_REAL;
  mean *= scale; second *= scale*scale;
  mse = second - 2*theta % mean + arma::square(theta);
  return Rcpp::List::create(Rcpp::Named("mean")=vec_r(mean), Rcpp::Named("second")=vec_r(second),
    Rcpp::Named("mse")=vec_r(mse), Rcpp::Named("risk")=arma::accu(mse),
    Rcpp::Named("theta")=vec_r(theta), Rcpp::Named("m")=vec_r(m), Rcpp::Named("s")=vec_r(s),
    Rcpp::Named("scale_opt")=copt,Rcpp::Named("risk_opt")=ropt);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_sure_grid(const arma::mat& x, const arma::vec& y, const arma::mat& b,
                         const arma::mat& v, const arma::vec& d, const arma::vec& k,
                         const arma::mat& gamma, double sigma2, bool intercept, double scale = 1.0) {
  if (x.n_rows != y.n_elem || x.n_cols != b.n_rows || b.n_cols != k.n_elem ||
      v.n_rows != b.n_rows || v.n_cols != d.n_elem || gamma.n_rows != k.n_elem)
    Rcpp::stop("incompatible SURE grid dimensions");
  if (!gamma.is_finite() || arma::any(arma::vectorise(gamma) < 0))
    Rcpp::stop("gamma must be finite and nonnegative");
  arma::mat sure(k.n_elem, gamma.n_cols), df(k.n_elem, gamma.n_cols), scale_min(k.n_elem,gamma.n_cols);
  arma::mat v2 = arma::square(v);
  arma::vec d2 = arma::square(d);
  for (arma::uword i = 0; i < k.n_elem; ++i) {
    Rcpp::checkUserInterrupt();
    arma::vec ak = v2 * (d2 / (d2 + k[i]));
    arma::mat th(b.n_rows, gamma.n_cols);
    for (arma::uword g = 0; g < gamma.n_cols; ++g) th.col(g) = threshold(b.col(i), gamma(i, g));
    arma::mat fits = x * th;
    arma::mat residual = scale * fits;
    residual.each_col() -= y;
    for (arma::uword g = 0; g < gamma.n_cols; ++g) {
      double df0 = 0;
      for (arma::uword j = 0; j < b.n_rows; ++j) if (th(j, g) != 0) df0 += ak[j];
      df(i,g) = scale*df0 + intercept;
      double den = arma::dot(fits.col(g),fits.col(g));
      scale_min(i,g) = den > 0 ? (arma::dot(y,fits.col(g))-sigma2*df0)/den : NA_REAL;
      sure(i, g) = arma::dot(residual.col(g), residual.col(g)) - x.n_rows * sigma2 + 2 * sigma2 * df(i, g);
    }
  }
  return Rcpp::List::create(Rcpp::Named("sure")=sure, Rcpp::Named("df")=df,Rcpp::Named("scale_min")=scale_min);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_mp_refit(const arma::mat& x, const arma::vec& y) {
  if (x.n_rows != y.n_elem) Rcpp::stop("incompatible refit dimensions");
  unsigned int rank;
  arma::vec b = refit(x, y, rank);
  return Rcpp::List::create(Rcpp::Named("coef")=vec_r(b), Rcpp::Named("rank")=static_cast<int>(rank),
                            Rcpp::Named("size")=static_cast<int>(x.n_cols));
}

// [[Rcpp::export(rng = false)]]
arma::mat cpp_cv_scores(const arma::mat& x, const arma::vec& y, const arma::mat& xtest,
                        const arma::vec& ytest, const arma::mat& b, const arma::mat& gamma,
                        int max_size, bool scale = false, bool hard = false) {
  if(max_size && (scale || hard)) Rcpp::stop("refit cannot be combined with scaling or hard thresholding");
  if (x.n_rows != y.n_elem || xtest.n_rows != ytest.n_elem ||
      x.n_cols != b.n_rows || xtest.n_cols != b.n_rows || max_size < 0 ||
      static_cast<arma::uword>(max_size) > b.n_rows)
    Rcpp::stop("incompatible CV dimensions");
  if (!max_size && (gamma.n_rows != b.n_cols || !gamma.is_finite() ||
      arma::any(arma::vectorise(gamma) < 0))) Rcpp::stop("invalid nonnegative threshold grid");
  arma::uword count = max_size > 0 ? static_cast<arma::uword>(max_size) : gamma.n_cols;
  arma::mat err(b.n_cols, count, arma::fill::zeros);
  for (arma::uword i = 0; i < b.n_cols; ++i) {
    Rcpp::checkUserInterrupt();
    if (max_size > 0) {
      std::vector<arma::uword> order(b.n_rows);
      std::iota(order.begin(), order.end(), 0);
      std::stable_sort(order.begin(), order.end(), [&](arma::uword a, arma::uword c) {
        return std::abs(b(a, i)) > std::abs(b(c, i));
      });
      for (arma::uword m = 1; m <= count; ++m) {
        arma::uvec selected(m);
        for (arma::uword j = 0; j < m; ++j) selected[j] = order[j];
        unsigned int rank;
        arma::vec cf = refit(x.cols(selected), y, rank);
        arma::vec residual = ytest - xtest.cols(selected) * cf;
        err(i, m - 1) = arma::dot(residual, residual);
      }
    } else {
      arma::mat th(b.n_rows, count);
      for (arma::uword g = 0; g < count; ++g) {
        th.col(g) = threshold_rule(b.col(i), gamma(i,g),hard);
        if(scale) th.col(g) *= ls_scale(x,y,th.col(g));
      }
      arma::mat residual = xtest * th;
      residual.each_col() -= ytest;
      err.row(i) = arma::sum(arma::square(residual), 0);
    }
  }
  return err;
}
