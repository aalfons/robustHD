\name{tslarsP}
\alias{print.tslarsP}
\alias{tslarsP}
\alias{tslarsP.default}
\alias{tslarsP.formula}
\title{Least angle regression for time series data with fixed lag length}
\usage{
  tslarsP(x, ...)

  \method{tslarsP}{formula} (formula, data, ...)

  \method{tslarsP}{default} (x, y, h = 1, p = 2, sMax = NA,
    fit = TRUE, crit = "BIC", ncores = 1, model = TRUE,
    ...)
}
\arguments{
  \item{formula}{a formula describing the full model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{tslarsP} is called.}

  \item{x}{a numeric matrix or data frame containing the
  candidate predictor series.}

  \item{y}{a numeric vector containing the response
  series.}

  \item{h}{an integer giving the forecast horizon (defaults
  to 1).}

  \item{p}{an integer giving the number of lags in the
  model (defaults to 2).}

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

  \item{fit}{a logical indicating whether to fit submodels
  along the sequence (\code{TRUE}, the default) or to
  simply return the sequence (\code{FALSE}).}

  \item{crit}{a character string specifying the optimality
  criterion to be used for selecting the final model.
  Currently, only \code{"BIC"} for the Bayes information
  criterion is implemented.}

  \item{ncores}{a positive integer giving the number of
  processor cores to be used for parallel computing (the
  default is 1 for no parallelization).  If this is set to
  \code{NA}, all available processor cores are used.}

  \item{model}{a logical indicating whether the model data
  should be included in the returned object.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  If \code{fit} is \code{FALSE}, an integer vector
  containing the indices of the sequenced predictor groups.

  Otherwise an object of class \code{"tslarsP"} (inheriting
  from class \code{"grplars"}) with the following
  components:

  \item{active}{an integer vector containing the sequence
  of predictor groups.}

  \item{coefficients}{a numeric matrix in which each column
  contains the regression coefficients of the corresponding
  submodel along the sequence.}

  \item{fitted.values}{a numeric matrix in which each
  column contains the fitted values of the corresponding
  submodel along the sequence.}

  \item{residuals}{a numeric matrix in which each column
  contains the residuals of the corresponding submodel
  along the sequence.}

  \item{crit}{a character string specifying the optimality
  criterion used for selecting the final model.}

  \item{critValues}{a numeric vector containing the values
  of the optimality criterion from the submodels along the
  sequence.}

  \item{df}{an integer vector containing the degrees of
  freedom of the submodels along the sequence (i.e., the
  number of estimated coefficients).}

  \item{sOpt}{an integer giving the optimal submodel.}

  \item{muY}{numeric; the mean of the response.}

  \item{sigmaY}{numeric; the standard deviation of the
  response.}

  \item{muX}{a numeric vector containing the means of the
  predictor variables.}

  \item{sigmaX}{a numeric vector containing the standard
  deviations of the predictor variables.}

  \item{x}{the matrix of candidate predictor series (if
  \code{model} is \code{TRUE}).}

  \item{y}{the response series (if \code{model} is
  \code{TRUE}).}

  \item{assign}{an integer vector giving the predictor
  group to which each predictor variable belongs.}

  \item{robust}{a logical indicating whether a robust fit
  was computed (\code{FALSE} for function \code{tslarsP}).}

  \item{h}{the forecast horizon.}

  \item{p}{the number of lags in the model.}

  \item{call}{the matched function call.}
}
\description{
  Sequence groups of candidate predictors and their
  respective lagged values according to their predictive
  content and find the optimal model along the sequence.
  Note that lagged values of the response are included as a
  predictor group as well.
}
\note{
  The predictor group of lagged values of the response is
  indicated by the index 0.

  Parallel computing for some of of the more
  computer-intensive computations for the inactive
  predictor groups is implemented on the C++ level via
  OpenMP (\url{http://openmp.org/}).
}
\author{
  Andreas Alfons, based on code by Sarah Gelper
}
\references{
  Gelper, S. and Croux, C. (2010) Time series least angle
  regression for selecting predictive economic sentiment
  series. Working paper.
}
\seealso{
  \code{\link{coef.seqModel}},
  \code{\link{fitted.seqModel}},
  \code{\link{residuals.seqModel}},
  \code{\link{predict.tslarsP}},
  \code{\link{plot.seqModel}}, \code{\link{tslars}},
  \code{\link{rtslarsP}}
}
\keyword{regression}
\keyword{ts}

