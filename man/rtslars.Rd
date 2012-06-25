\name{rtslars}
\alias{rtslars}
\alias{rtslars.default}
\alias{rtslars.formula}
\title{Robust least angle regression for time series data}
\usage{
  rtslars(x, ...)

  \method{rtslars}{formula} (formula, data, ...)

  \method{rtslars}{default} (x, y, h = 1, pMax = 3,
    sMax = NA, centerFun = median, scaleFun = mad,
    regFun = lmrob, regArgs = list(), winsorize = FALSE,
    const = 2, prob = 0.95,
    combine = c("min", "euclidean", "mahalanobis"),
    ncores = 1, cl = NULL, fit = TRUE, crit = "BIC",
    model = TRUE, ...)
}
\arguments{
  \item{formula}{a formula describing the full model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{rtslars} is called.}

  \item{x}{a numeric matrix or data frame containing the
  candidate predictor series.}

  \item{y}{a numeric vector containing the response
  series.}

  \item{h}{an integer giving the forecast horizon (defaults
  to 1).}

  \item{pMax}{an integer giving the maximum number of lags
  in the model (defaults to 3).}

  \item{sMax}{an integer vector of length two.  If a single
  integer is supplied, it is recycled.  The first element
  gives the number of predictor groups to be sequenced.  If
  it is \code{NA} (the default), predictor groups are
  sequenced as long as there are no singularity issues.
  The second element gives the maximum number of predictor
  groups to be included in the final model.  If it is
  \code{NA} (the default), predictor groups may be added to
  the model as long as there are twice as many observations
  as predictor variables.}

  \item{centerFun}{a function to compute a robust estimate
  for the center (defaults to
  \code{\link[stats]{median}}).}

  \item{scaleFun}{a function to compute a robust estimate
  for the scale (defaults to \code{\link[stats]{mad}}).}

  \item{regFun}{a function to compute robust linear
  regressions that can be interpreted as weighted least
  squares (defaults to \code{\link[robustbase]{lmrob}}).}

  \item{regArgs}{a list of arguments to be passed to
  \code{regFun}.}

  \item{winsorize}{a logical indicating whether to clean
  the data by multivariate winsorization.}

  \item{const}{numeric; tuning constant for multivariate
  winsorization to be used in the initial corralation
  estimates based on adjusted univariate winsorization
  (defaults to 2).}

  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution to be used in
  multivariate winsorization (defaults to 0.95).}

  \item{combine}{a character string specifying how to
  combine the data cleaning weights from the robust
  regressions with each predictor group. Possible values
  are \code{"min"} for taking the minimum weight for each
  observation, \code{"euclidean"} for weights based on
  Euclidean distances of the multivariate set of
  standardized residuals (i.e., multivariate winsorization
  of the standardized residuals assuming independence), or
  \code{"mahalanobis"} for weights based on Mahalanobis
  distances of the multivariate set of standardized
  residuals (i.e., multivariate winsorization of the
  standardized residuals).}

  \item{fit}{a logical indicating whether to fit submodels
  along the sequence (\code{TRUE}, the default) or to
  simply return the sequence (\code{FALSE}).}

  \item{crit}{a character string specifying the optimality
  criterion to be used for selecting the final model.
  Currently, only \code{"BIC"} for the Bayes information
  criterion based on a robust residual scale estimate is
  implemented.}

  \item{model}{a logical indicating whether the model data
  should be included in the returned object.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  If \code{fit} is \code{FALSE}, an integer matrix in which
  each column contains the indices of the sequenced
  predictor groups for the corresponding lag length.

  Otherwise an object of class \code{"tslars"} with the
  following components:

  \item{pFit}{a list containing the fits for the respective
  lag lengths. (see \code{\link{rtslarsP}}).}

  \item{pOpt}{an integer giving the optimal number of
  lags.}

  \item{pMax}{the maximum number of lags considered.}

  \item{x}{the matrix of candidate predictor series (if
  \code{model} is \code{TRUE}).}

  \item{y}{the response series (if \code{model} is
  \code{TRUE}).}

  \item{call}{the matched function call.}
}
\description{
  Robustly sequence groups of candidate predictors and
  their respective lagged values according to their
  predictive content and find the optimal model along the
  sequence.  Note that lagged values of the response are
  included as a predictor group as well.
}
\note{
  The predictor group of lagged values of the response is
  indicated by the index 0.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{coef.tslars}}, \code{\link{fitted.tslars}},
  \code{\link{residuals.tslars}},
  \code{\link{predict.tslars}}, \code{\link{plot.tslars}},
  \code{\link{rtslarsP}}, \code{\link{tslars}}
}
\keyword{regression}
\keyword{robust}
\keyword{ts}

