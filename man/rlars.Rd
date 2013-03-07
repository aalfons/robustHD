\name{rlars}
\alias{print.rlars}
\alias{rlars}
\alias{rlars.default}
\alias{rlars.formula}
\title{Robust least angle regression}
\usage{
  rlars(x, ...)

  \method{rlars}{formula} (formula, data, ...)

  \method{rlars}{default} (x, y, sMax = NA,
    centerFun = median, scaleFun = mad, winsorize = FALSE,
    pca = FALSE, const = 2, prob = 0.95, fit = TRUE,
    s = c(0, sMax), regFun = lmrob, regArgs = list(),
    crit = c("BIC", "PE"), splits = foldControl(),
    cost = rtmspe, costArgs = list(),
    selectBest = c("hastie", "min"), seFactor = 1,
    ncores = 1, cl = NULL, seed = NULL, model = TRUE,
    tol = .Machine$double.eps^0.5, ...)
}
\arguments{
  \item{formula}{a formula describing the full model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{rlars} is called.}

  \item{x}{a matrix or data frame containing the candidate
  predictors.}

  \item{y}{a numeric vector containing the response.}

  \item{sMax}{an integer giving the number of predictors to
  be sequenced.  If it is \code{NA} (the default),
  predictors are sequenced as long as there are twice as
  many observations as predictors.}

  \item{centerFun}{a function to compute a robust estimate
  for the center (defaults to
  \code{\link[stats]{median}}).}

  \item{scaleFun}{a function to compute a robust estimate
  for the scale (defaults to \code{\link[stats]{mad}}).}

  \item{winsorize}{a logical indicating whether to clean
  the full data set by multivariate winsorization, i.e., to
  perform data cleaning RLARS instead of plug-in RLARS
  (defaults to \code{FALSE}).}

  \item{pca}{a logical indicating whether a robust PCA step
  should be performed when computing the data cleaning
  weights for multivariate winsorization (defaults to
  \code{FALSE}).  The distances of the observations are
  then computed on the PCA scores rather than the original
  observations, which makes data cleaning RLARS applicable
  for high-dimensional data.  If \code{TRUE} or \code{NA},
  components are computed as long as the robust correlation
  matrix of the scores can be inverted, and all components
  with an eigenvalue larger than or equal to 1 are
  retained.  Alternatively, an integer giving the maximum
  number of components can be supplied.}

  \item{const}{numeric; tuning constant to be used in the
  initial corralation estimates based on adjusted
  univariate winsorization (defaults to 2).}

  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution to be used in
  bivariate or multivariate winsorization (defaults to
  0.95).}

  \item{fit}{a logical indicating whether to fit submodels
  along the sequence (\code{TRUE}, the default) or to
  simply return the sequence (\code{FALSE}).}

  \item{s}{an integer vector of length two giving the first
  and last step along the sequence for which to compute
  submodels.  The default is to start with a model
  containing only an intercept (step 0) and iteratively add
  all variables along the sequence (step \code{sMax}).  If
  the second element is \code{NA}, predictors may be added
  to the model as long as there are twice as many
  observations as predictors.  If only one value is
  supplied, it is used as the last step, with the first
  step being 0.}

  \item{regFun}{a function to compute robust linear
  regressions along the sequence (defaults to
  \code{\link[robustbase]{lmrob}}).}

  \item{regArgs}{a list of arguments to be passed to
  \code{regFun}.}

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
  fitting models along the sequence and for prediction
  error estimation, parallel computing is implemented on
  the \R level using package \pkg{parallel}.  Otherwise
  parallel computing for some of of the more
  computer-intensive computations in the sequencing step is
  implemented on the C++ level via OpenMP
  (\url{http://openmp.org/}).}

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

  \item{tol}{a small positive numeric value.  This is used
  in bivariate winsorization to determine whether the
  initial estimate from adjusted univariate winsorization
  is close to 1 in absolute value.  In this case, bivariate
  winsorization would fail since the points form almost a
  straight line, and the initial estimate is returned.}

  \item{\dots}{additional arguments to be passed down.  For
  the default method, additional arguments to be passed
  down to \code{\link[=standardize]{robStandardize}}.}
}
\value{
  If \code{fit} is \code{FALSE}, an integer vector
  containing the indices of the sequenced predictors.

  Otherwise an object of class \code{"rlars"} (inheriting
  from class \code{"seqModel"}) with the following
  components:

  \item{active}{an integer vector containing the indices of
  the sequenced predictors.}

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
  was computed (\code{TRUE} for \code{"rlars"} models).}

  \item{scale}{a numeric vector giving the robust residual
  scale estimates for the submodels along the sequence.}

  \item{crit}{an object of class \code{"bicSelect"}
  containing the BIC values and indicating the final model
  (only returned if argument \code{crit} is \code{"BIC"}
  and argument \code{s} indicates more than one step along
  the sequence).}

  \item{muX}{a numeric vector containing the center
  estimates of the predictors.}

  \item{sigmaX}{a numeric vector containing the scale
  estimates of the predictors.}

  \item{muY}{numeric; the center estimate of the response.}

  \item{sigmaY}{numeric; the scale estimate of the
  response.}

  \item{x}{the matrix of candidate predictors (if
  \code{model} is \code{TRUE}).}

  \item{y}{the response (if \code{model} is \code{TRUE}).}

  \item{w}{a numeric vector giving the data cleaning
  weights (if \code{winsorize} is \code{TRUE}).}

  \item{call}{the matched function call.}
}
\description{
  Robustly sequence candidate predictors according to their
  predictive content and find the optimal model along the
  sequence.
}
\examples{
## generate data
# example is not high-dimensional to keep computation time low
library("mvtnorm")
set.seed(1234)  # for reproducibility
n <- 100  # number of observations
p <- 25   # number of variables
beta <- rep.int(c(1, 0), c(5, p-5))  # coefficients
sigma <- 0.5      # controls signal-to-noise ratio
epsilon <- 0.1    # contamination level
Sigma <- 0.5^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
x <- rmvnorm(n, sigma=Sigma)    # predictor matrix
e <- rnorm(n)                   # error terms
i <- 1:ceiling(epsilon*n)       # observations to be contaminated
e[i] <- e[i] + 5                # vertical outliers
y <- c(x \%*\% beta + sigma * e)  # response
x[i,] <- x[i,] + 5              # bad leverage points

## fit robust LARS model
rlars(x, y)
}
\author{
  Andreas Alfons, based on code by Jafar A. Khan, Stefan
  Van Aelst and Ruben H. Zamar
}
\references{
  Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust
  linear model selection based on least angle regression.
  \emph{Journal of the American Statistical Association},
  \bold{102}(480), 1289--1299.
}
\seealso{
  \code{\link{coef.seqModel}},
  \code{\link{fitted.seqModel}},
  \code{\link{residuals.seqModel}},
  \code{\link{predict.seqModel}},
  \code{\link{plot.seqModel}}
}
\keyword{regression}
\keyword{robust}

