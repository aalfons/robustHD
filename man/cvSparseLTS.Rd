\name{cvSparseLTS}
\alias{cvSparseLTS}
\alias{cvSparseLTS.default}
\alias{cvSparseLTS.formula}
\alias{print.cvSparseLTS}
\title{Cross-validation for sparse LTS regression models}
\usage{
  cvSparseLTS(x, ...)

  \method{cvSparseLTS}{formula} (formula, data, ...)

  \method{cvSparseLTS}{default} (x, y, cost = rtmspe,
    K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, fit = c("reweighted", "raw", "both"),
    selectBest = c("min", "hastie"), seFactor = 1, lambda,
    mode = c("lambda", "fraction"), alpha = 0.75,
    intercept = TRUE, nsamp = c(500, 10),
    initial = c("sparse", "hyperplane", "random"),
    ncstep = 2, use.correction = TRUE,
    tol = .Machine$double.eps^0.5,
    eps = .Machine$double.eps, use.Gram = TRUE,
    centerFun = median, scaleFun = mad, const = 2,
    prob = 0.95, fallback = FALSE, seed = NULL, ...)
}
\arguments{
  \item{formula}{a formula describing the model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{cvSparseLTS} is called.}

  \item{x}{a numeric matrix containing the predictor
  variables.}

  \item{y}{a numeric vector containing the response
  variable.}

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

  \item{selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter selects the model with the largest
  penalty parameter whose prediction error is no larger
  than \code{seFactor} standard errors above the prediction
  error of the best overall model.}

  \item{seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{selectBest} is
  \code{"min"}.}

  \item{lambda}{a numeric vector of non-negative values
  giving the penalty parameters for which to estimate the
  prediction error.}

  \item{mode}{a character string specifying the type of
  penalty parameter.  If \code{"lambda"}, this gives
  penalty parameter directly.  If \code{"fraction"}, the
  smallest value of the penalty parameter that sets all
  coefficients to 0 is first estimated based on bivariate
  winsorization, then \code{lambda} gives the fraction of
  that estimate to be used (hence \code{lambda} should be
  in the interval [0,1] in that case).}

  \item{alpha}{a numeric value giving the percentage of the
  residuals for which the \eqn{L_{1}}{L1} penalized sum of
  squares should be minimized (the default is 0.75).}

  \item{intercept}{a logical indicating whether a constant
  term should be included in the model (the default is
  \code{TRUE}).}

  \item{nsamp}{a numeric vector giving the number of
  subsamples to be used in the two phases of the algorithm.
  The first element gives the number of initial subsamples
  to be used.  The second element gives the number of
  subsamples to keep after the first phase of \code{ncstep}
  C-steps.  For those remaining subsets, additional C-steps
  are performed until convergence.  The default is to first
  perform \code{ncstep} C-steps on 500 initial subsamples,
  and then to keep the 10 subsamples with the lowest value
  of the objective function for additional C-steps until
  convergence.}

  \item{initial}{a character string specifying the type of
  initial subsamples to be used.  If \code{"sparse"}, the
  lasso fit given by three randomly selected data points is
  first computed.  The corresponding initial subsample is
  then formed by the fraction \code{alpha} of data points
  with the smallest squared residuals.  Note that this is
  optimal from a robustness point of view, as the
  probability of including an outlier in the initial lasso
  fit is minimized.  If \code{"hyperplane"}, a hyperplane
  through \eqn{p} randomly selected data points is first
  computed, where \eqn{p} denotes the number of variables.
  The corresponding initial subsample is then again formed
  by the fraction \code{alpha} of data points with the
  smallest squared residuals. Note that this cannot be
  applied if \eqn{p} is larger than the number of
  observations.  Nevertheless, the probability of including
  an outlier increases with increasing dimension \eqn{p}.
  If \code{"random"}, the initial subsamples are given by a
  fraction \code{alpha} of randomly selected data points.
  Note that this leads to the largest probability of
  including an outlier.}

  \item{ncstep}{a positive integer giving the number of
  C-steps to perform on all subsamples in the first phase
  of the algorithm (the default is to perform two
  C-steps).}

  \item{use.correction}{currently ignored.  The consistency
  factor for the residual scale estimate is always
  applied.}

  \item{tol}{a small positive numeric value giving the
  tolerance for convergence.}

  \item{eps}{a small positive numeric value used to
  determine whether the variability within a variable is
  too small (an effective zero).}

  \item{use.Gram}{a logical indicating whether the Gram
  matrix of the explanatory variables should be precomputed
  in the lasso fits (the default is \code{TRUE}).  If the
  number of variables is large (e.g., larger than the
  number of observations), computation may be faster when
  this is set to \code{FALSE}.}

  \item{centerFun}{a function to compute a robust estimate
  for the center to be used for robust standardization
  (defaults to \code{\link[stats]{median}}).  Ignored if
  \code{standardized} is \code{TRUE}.}

  \item{scaleFun}{a function to compute a robust estimate
  for the scale to be used for robust standardization
  (defaults to \code{\link[stats]{mad}}). Ignored if
  \code{standardized} is \code{TRUE}.}

  \item{const}{numeric; tuning constant to be used in
  univariate winsorization (defaults to 2).}

  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution to be used in
  multivariate winsorization (defaults to 0.95).}

  \item{fallback}{a logical indicating whether
  standardization with \code{\link{mean}} and
  \code{\link[stats]{sd}} should be performed as a fallback
  mode for variables whose robust scale estimate is too
  small.  This is useful, e.g., for data containing dummy
  variables.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed to the
  prediction loss function \code{cost}.}
}
\value{
  An object of class \code{"cvSparseLTS"} (which inherits
  from classes \code{"cvTuning"} and \code{"cvSelect"})
  with the following components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{tuning}{a data frame containing the values of the
  penalty parameter for which the prediction error was
  estimated.}

  \item{best}{an integer vector giving the indices of the
  optimal penalty parameter for the requested model fits.}

  \item{cv}{a data frame containing the estimated
  prediction errors of the requested model fits for all
  values of the penalty parameter.}

  \item{se}{a data frame vector containing the estimated
  standard errors of the prediction loss for the requested
  model fits for all values of the penalty parameter.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{reps}{a data frame containing the estimated
  prediction errors of the requested model fits from all
  replications and for all values of the penalty parameter.
  This is only returned for repeated cross-validation.}

  \item{seed}{the seed of the random number generator
  before cross-validation was performed.}

  \item{call}{the matched function call.}
}
\description{
  Estimate the prediction error of sparse least trimmed
  squares regression models for various values of the
  penalty parameter via (repeated) \eqn{K}-fold
  cross-validation.
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

## evaluate sparse LTS models over a grid of values for lambda
frac <- seq(0.25, 0.05, by = -0.05)
cv <- cvSparseLTS(x, y, lambda = frac, mode = "fraction", 
    selectBest = "hastie", includeSE = TRUE)
cv
plot(cv)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{repCV.sparseLTS}}, \code{\link{sparseLTS}},
  \code{\link{predict.sparseLTS}},
  \code{\link[cvTools]{cvFolds}},
  \code{\link[cvTools]{cost}}
}
\keyword{utilities}

