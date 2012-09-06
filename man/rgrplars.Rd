\name{rgrplars}
\alias{rgrplars}
\alias{rgrplars.data.frame}
\alias{rgrplars.default}
\alias{rgrplars.formula}
\title{Robust groupwise least angle regression}
\usage{
  rgrplars(x, ...)

  \method{rgrplars}{formula} (formula, data, ...)

  \method{rgrplars}{data.frame} (x, y, ...)

  \method{rgrplars}{default} (x, y, sMax = NA, assign,
    dummy, centerFun = median, scaleFun = mad,
    regFun = lmrob, regArgs = list(),
    combine = c("min", "euclidean", "mahalanobis"),
    const = 2, prob = 0.95, fit = TRUE,
    crit = c("BIC", "PE"), splits = foldControl(),
    cost = rtmspe, costArgs = list(),
    selectBest = c("hastie", "min"), seFactor = 1,
    ncores = 1, cl = NULL, seed = NULL, model = TRUE, ...)
}
\arguments{
  \item{formula}{a formula describing the full model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{rgrplars} is called.}

  \item{x}{a matrix or data frame containing the candidate
  predictors.}

  \item{y}{a numeric vector containing the response.}

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

  \item{assign}{an integer vector giving the predictor
  group to which each predictor variable belongs.}

  \item{dummy}{a logical vector indicating whether the
  predictors are dummy variables.}

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

  \item{const}{numeric; tuning constant for multivariate
  winsorization to be used in the initial corralation
  estimates based on adjusted univariate winsorization
  (defaults to 2).}

  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution to be used in
  multivariate winsorization (defaults to 0.95).}

  \item{fit}{a logical indicating whether to fit submodels
  along the sequence (\code{TRUE}, the default) or to
  simply return the sequence (\code{FALSE}).}

  \item{crit}{a character string specifying the optimality
  criterion to be used for selecting the final model.
  Possible values are \code{"BIC"} for the Bayes
  information criterion and \code{"PE"} for
  resampling-based prediction error estimation.}

  \item{splits}{an object giving data splits to be used for
  prediction error estimation (see
  \code{\link[perry]{perry}}).}

  \item{cost}{a cost function measuring prediction loss
  (see \code{\link[perry]{perry}} for some requirements).
  The default is to use the root trimmed mean squared
  prediction error (see \code{\link[perry]{cost}}).}

  \item{costArgs}{a list of additional arguments to be
  passed to the prediction loss function \code{cost}.}

  \item{selectBest,seFactor}{arguments specifying a
  criterion for selecting the best model (see
  \code{\link[perry]{perrySelect}}).  The default is to use
  a one-standard-error rule.}

  \item{ncores}{a positive integer giving the number of
  processor cores to be used for parallel computing (the
  default is 1 for no parallelization).  If this is set to
  \code{NA}, all available processor cores are used.  For
  obtaining the data cleaning weights, for fitting models
  along the sequence and for prediction error estimation,
  parallel computing is implemented on the \R level using
  package \pkg{parallel}.  Otherwise parallel computing for
  some of of the more computer-intensive computations in
  the sequencing step is implemented on the C++ level via
  OpenMP (\url{http://openmp.org/}).}

  \item{cl}{a \pkg{parallel} cluster for parallel computing
  as generated by \code{\link[parallel]{makeCluster}}.
  This is preferred over \code{ncores} for tasks that are
  parallelized on the \R level, in which case \code{ncores}
  is only used for tasks that are parallelized on the C++
  level.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).  This is
  useful because many robust regression functions
  (including \code{\link[robustbase]{lmrob}}) involve
  randomness, or for prediction error estimation.  On
  parallel \R worker processes, random number streams are
  used and the seed is set via
  \code{\link{clusterSetRNGStream}}.}

  \item{model}{a logical indicating whether the model data
  should be included in the returned object.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  If \code{fit} is \code{FALSE}, an integer vector
  containing the indices of the sequenced predictor groups.

  Otherwise an object of class \code{"grplars"} (inheriting
  from class \code{"seqModel"} if \code{crit="BIC"} or
  \code{"optSeqModel"} if \code{crit="PE"}) with the
  following components:

  \item{active}{an integer vector containing the sequence
  of predictor groups.}

  \item{df}{an integer vector containing the degrees of
  freedom of the submodels along the sequence (i.e., the
  number of estimated coefficients).}

  \item{coefficients}{a numeric matrix in which each column
  contains the regression coefficients of the corresponding
  submodel along the sequence (\code{"seqModel"}); or a
  numeric vector containing the coefficients of the optimal
  submodel (\code{"optSeqModel"}).}

  \item{fitted.values}{a numeric matrix in which each
  column contains the fitted values of the corresponding
  submodel along the sequence (\code{"seqModel"}); or a
  numeric vector containing the fitted values of the
  optimal submodel (\code{"optSeqModel"}).}

  \item{residuals}{a numeric matrix in which each column
  contains the residuals of the corresponding submodel
  along the sequence (\code{"seqModel"}); or a numeric
  vector containing the residuals of the optimal submodel
  (\code{"optSeqModel"}).}

  \item{crit}{a character string specifying the optimality
  criterion used for selecting the final model.}

  \item{critValues}{a numeric vector containing the values
  of the optimality criterion from the submodels along the
  sequence (\code{"seqModel"}); or an object of class
  \code{"perrySeqModel"} (inheriting from
  \code{"\link[perry]{perrySelect}"}) that contains the
  estimated prediction errors of the submodels
  (\code{"optSeqModel"}).}

  \item{sOpt}{an integer giving the optimal submodel (only
  \code{"seqModel"}).}

  \item{muY}{numeric; the mean of the cleaned response.}

  \item{sigmaY}{numeric; the standard deviation of the
  cleaned response.}

  \item{muX}{a numeric vector containing the means of the
  cleaned predictor variables.}

  \item{sigmaX}{a numeric vector containing the standard
  deviations of the cleaned predictor variables.}

  \item{x}{the matrix of candidate predictors (if
  \code{model} is \code{TRUE}).}

  \item{y}{the response (if \code{model} is \code{TRUE}).}

  \item{assign}{an integer vector giving the predictor
  group to which each predictor variable belongs.}

  \item{robust}{a logical indicating whether a robust fit
  was computed (\code{TRUE} for function \code{rgrplars}).}

  \item{w}{a numeric vector giving the data cleaning
  weights.}

  \item{call}{the matched function call.}
}
\description{
  Robustly sequence groups of candidate predictors
  according to their predictive content and find the
  optimal model along the sequence.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{coef.seqModel}},
  \code{\link{fitted.seqModel}},
  \code{\link{residuals.seqModel}},
  \code{\link{predict.seqModel}},
  \code{\link{plot.seqModel}}, \code{\link{grplars}}
}
\keyword{regression}

