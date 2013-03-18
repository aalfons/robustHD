\name{grplars}
\alias{grplars}
\alias{grplars.data.frame}
\alias{grplars.default}
\alias{grplars.formula}
\alias{print.grplars}
\alias{rgrplars}
\alias{rgrplars.data.frame}
\alias{rgrplars.default}
\alias{rgrplars.formula}
\title{(Robust) groupwise least angle regression}
\usage{
  grplars(x, ...)

  \method{grplars}{formula} (formula, data, ...)

  \method{grplars}{data.frame} (x, y, ...)

  \method{grplars}{default} (x, y, sMax = NA, assign,
    fit = TRUE, s = c(0, sMax), crit = c("BIC", "PE"),
    splits = foldControl(), cost = rmspe,
    costArgs = list(), selectBest = c("hastie", "min"),
    seFactor = 1, ncores = 1, cl = NULL, seed = NULL,
    model = TRUE, ...)

  rgrplars(x, ...)

  \method{rgrplars}{formula} (formula, data, ...)

  \method{rgrplars}{data.frame} (x, y, ...)

  \method{rgrplars}{default} (x, y, sMax = NA, assign,
    dummy, centerFun = median, scaleFun = mad,
    regFun = lmrob, regArgs = list(),
    combine = c("min", "mahalanobis"), pca = FALSE,
    const = 2, prob = 0.95, fit = TRUE, s = c(0, sMax),
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
  environment from which \code{grplars} or \code{rgrplars}
  is called.}

  \item{x}{a matrix or data frame containing the candidate
  predictors.}

  \item{y}{a numeric vector containing the response.}

  \item{sMax}{an integer giving the number of predictor
  groups to be sequenced.  If it is \code{NA} (the
  default), predictor groups are sequenced as long as there
  are twice as many observations as predictor variables.}

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
  observation, or \code{"mahalanobis"} for weights based on
  Mahalanobis distances of the multivariate set of
  standardized residuals (i.e., multivariate winsorization
  of the standardized residuals).}

  \item{pca}{a logical indicating whether a robust PCA step
  should be performed when computing the data cleaning
  weights based on Mahalanobis distances of the
  standardized residuals (defaults to \code{FALSE}).  The
  distances of the observations are then computed on the
  PCA scores, which makes the approach applicable in the
  high-dimensional case.  If \code{TRUE} or \code{NA},
  components are computed as long as the robust correlation
  matrix of the scores can be inverted, and all components
  with an eigenvalue larger than or equal to 1 are
  retained.  Alternatively, an integer giving the maximum
  number of components can be supplied.}

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

  \item{s}{an integer vector of length two giving the first
  and last step along the sequence for which to compute
  submodels.  The default is to start with a model
  containing only an intercept (step 0) and iteratively add
  all groups along the sequence (step \code{sMax}).  If the
  second element is \code{NA}, predictor groups may be
  added to the model as long as there are twice as many
  observations as predictor variables.  If only one value
  is supplied, it is recycled.}

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
  prediction error for robust fits and the root mean
  squared prediction error otherwise (see
  \code{\link[perry]{cost}}).}

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

  Else if \code{crit} is \code{"PE"}, an object of class
  \code{"perrySeqModel"} (inheriting from classes
  \code{"perryTuning"}, see
  \code{\link[perry]{perryTuning}}).  It contains
  information on the prediction error criterion, and
  includes the final model as component \code{finalModel}.

  Otherwise an object of class \code{"grplars"} (inheriting
  from class \code{"seqModel"}) with the following
  components:

  \item{active}{an integer vector containing the sequence
  of predictor groups.}

  \item{s}{an integer vector containing the steps for which
  submodels along the sequence have been computed.}

  \item{coefficients}{a numeric matrix in which each column
  contains the regression coefficients of the corresponding
  submodel along the sequence.}

  \item{fitted.values}{a numeric matrix in which each
  column contains the fitted values of the corresponding
  submodel along the sequence.}

  \item{residuals}{a numeric matrix in which each column
  contains the residuals of the corresponding submodel
  along the sequence.}

  \item{df}{an integer vector containing the degrees of
  freedom of the submodels along the sequence (i.e., the
  number of estimated coefficients).}

  \item{robust}{a logical indicating whether a robust fit
  was computed.}

  \item{scale}{a numeric vector giving the robust residual
  scale estimates for the submodels along the sequence
  (only returned for robust fits).}

  \item{crit}{an object of class \code{"bicSelect"}
  containing the BIC values and indicating the final model
  (only returned if argument \code{crit} is \code{"BIC"}
  and argument \code{s} indicates more than one step along
  the sequence).}

  \item{muX}{a numeric vector containing the center
  estimates of the predictor variables.}

  \item{sigmaX}{a numeric vector containing the scale
  estimates of the predictor variables.}

  \item{muY}{numeric; the center estimate of the response.}

  \item{sigmaY}{numeric; the scale estimate of the
  response.}

  \item{x}{the matrix of candidate predictors (if
  \code{model} is \code{TRUE}).}

  \item{y}{the response (if \code{model} is \code{TRUE}).}

  \item{assign}{an integer vector giving the predictor
  group to which each predictor variable belongs.}

  \item{w}{a numeric vector giving the data cleaning
  weights (only returned for robust fits).}

  \item{call}{the matched function call.}
}
\description{
  (Robustly) sequence groups of candidate predictors
  according to their predictive content and find the
  optimal model along the sequence.
}
\note{
  \code{findStepSizes} is a utility function that computes
  the step size for each inactive predictor group.  It is
  only exported so it can be called by the underlying C++
  code for sequencing the predictor groups.  Hence it is
  not expected to be called by the user and not documented.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[=coef.seqModel]{coef}},
  \code{\link[=fitted.seqModel]{fitted}},
  \code{\link[=plot.seqModel]{plot}},
  \code{\link[=predict.seqModel]{predict}},
  \code{\link[=residuals.seqModel]{residuals}},
  \code{\link[robustbase]{lmrob}}
}
\keyword{regression}
\keyword{robust}

