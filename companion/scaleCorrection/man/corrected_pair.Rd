\name{corrected_pair}
\alias{corrected_pair}
\alias{estimate_variance}
\alias{variance_folds}
\title{Research implementation of a nonnegative optimism-corrected scale}
\description{Reuses frozen ridge preprocessing and ridge paths. For centered training
data, computes the scale as the positive part of the inner product of the
response with the unscaled fit minus estimated variance times the selected
ridge diagonal sum, divided by the squared norm of the unscaled fit.
Zero fitted vectors receive scale zero. This is a plug-in correction, not an
unbiased risk estimate for the adaptively selected estimator.}
\usage{
corrected_pair(X, y, grid, foldid, unit_id, variance_fun = estimate_variance)
estimate_variance(X, y, unit_id, subset = "full")
variance_folds(n, unit_id, subset = "full")
}
\arguments{
\item{X}{Finite numeric design matrix.}
\item{y}{Finite response vector.}
\item{grid}{List containing positive ridge penalties and threshold fractions, starting at zero.}
\item{foldid}{Saved outer cross-validation folds.}
\item{unit_id}{Frozen replication identifier.}
\item{variance_fun}{Training-subset variance estimator, replaceable for validation.}
\item{subset}{Deterministic namespace for full training or an outer training subset.}
\item{n}{Number of observations in a variance-estimation subset, at least 20.}
}
\value{The two selected fits, a shared CV surface and variance diagnostics; or
variance diagnostics or deterministic fold IDs for the respective helpers.}
\details{Variance uses ten-fold cross-validated Lasso at lambda.min, with
residual sum of squares divided by max(1,n-selected size-1), matching the
intercept adjustment in the supplied edition-4 code. It is recomputed within
each outer training subset. The original outer folds and data seeds are unchanged.
The fixed-threshold ridge diagonal expression is used as a plug-in term when
thresholds depend on the response. No new risk dominance claim is made.}
\references{Reid, Tibshirani and Friedman (2016), A Study of Error Variance
Estimation in Lasso Regression. Statistica Sinica 26, 35--67.}
