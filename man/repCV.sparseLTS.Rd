\name{repCV.sparseLTS}
\alias{repCV.sparseLTS}
\alias{repCV.sparseLTSGrid}
\title{Cross-validation for sparse LTS regression models}
\usage{
  \method{repCV}{sparseLTS} (object, cost = rtmspe, K = 5,
    R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, fit = c("reweighted", "raw", "both"),
    seed = NULL, ...)

  \method{repCV}{sparseLTSGrid} (object, cost = rtmspe,
    K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, fit = c("reweighted", "raw", "both"),
    seed = NULL, ...)
}
\arguments{
  \item{object}{the model fit for which to estimate the
  prediction error.}

  \item{cost}{a robust cost function measuring prediction
  loss.  It should expect vectors to be passed as its first
  two arguments, the first corresponding to the observed
  values of the response and the second to the predicted
  values, and must return a non-negative scalar value.  The
  default is to use the root trimmed mean squared
  prediction error (see \code{\link[cvTools]{cost}}).}

  \item{K}{an integer giving the number of groups into
  which the data should be split (the default is five).
  Keep in mind that this should be chosen such that all
  groups are of approximately equal size.  Setting \code{K}
  equal to \code{n} yields leave-one-out cross-validation.}

  \item{R}{an integer giving the number of replications for
  repeated \eqn{K}-fold cross-validation.  This is ignored
  for for leave-one-out cross-validation and other
  non-random splits of the data.}

  \item{foldType}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}

  \item{folds}{an object of class \code{"cvFolds"} giving
  the folds of the data for cross-validation (as returned
  by \code{\link[cvTools]{cvFolds}}).  If supplied, this is
  preferred over \code{K} and \code{R}.}

  \item{fit}{a character string specifying for which fit to
  estimate the prediction error.  Possible values are
  \code{"reweighted"} (the default) for the prediction
  error of the reweighted fit, \code{"raw"} for the
  prediction error of the raw fit, or \code{"both"} for the
  prediction error of both fits.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed to the
  prediction loss function \code{cost}.}
}
\value{
  An object of class \code{"cv"} with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{cv}{a numeric vector containing the estimated
  prediction errors for the requested model fits.  For
  repeated cross-validation, this contains the average
  values over all replications.}

  \item{se}{a numeric vector containing the estimated
  standard errors of the prediction loss for the requested
  model fits.}

  \item{reps}{a numeric matrix in which each column
  contains the estimated prediction errors from all
  replications for the requested model fits.  This is only
  returned for repeated cross-validation.}

  \item{seed}{the seed of the random number generator
  before cross-validation was performed.}

  \item{call}{the matched function call.}
}
\description{
  Estimate the prediction error of a previously fit sparse
  least trimmed squares regression model via (repeated)
  \eqn{K}-fold cross-validation.  If the model fit contains
  estimates over a grid of values for the penalty parameter
  (i.e., for objects of class \code{"sparseLTSGrid"}), then
  in each iteration of cross-validation the optimal model
  is selected from the training data and used to make
  predictions for the test data.
}
\examples{
## generate data
library("mvtnorm")
set.seed(1234)  # for reproducibility
n <- 100  # number of observations
p <- 200  # number of variables
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

## fit and evaluate sparse LTS model
fit <- sparseLTS(x, y, lambda = 0.05, mode = "fraction")
cv <- repCV(fit)
cv
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{sparseLTS}}, \code{\link{sparseLTSGrid}},
  \code{\link{predict.sparseLTS}},
  \code{\link[cvTools]{cvFolds}},
  \code{\link[cvTools]{cost}}
}
\keyword{robust}
\keyword{utilities}

