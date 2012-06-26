\name{sparseLTS}
\alias{print.sparseLTS}
\alias{sparseLTS}
\alias{sparseLTS.default}
\alias{sparseLTS.formula}
\title{Sparse least trimmed squares regression}
\usage{
  sparseLTS(x, ...)

  \method{sparseLTS}{formula} (formula, data, ...)

  \method{sparseLTS}{default} (x, y, lambda,
    mode = c("lambda", "fraction"), alpha = 0.75,
    intercept = TRUE, nsamp = c(500, 10),
    initial = c("sparse", "hyperplane", "random"),
    ncstep = 2, use.correction = TRUE,
    tol = .Machine$double.eps^0.5,
    eps = .Machine$double.eps, use.Gram = TRUE, ncores = 1,
    seed = NULL, model = TRUE, ...)
}
\arguments{
  \item{formula}{a formula describing the model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{sparseLTS} is called.}

  \item{x}{a numeric matrix containing the predictor
  variables.}

  \item{y}{a numeric vector containing the response
  variable.}

  \item{lambda}{a non-negative numeric value giving the
  penalty parameter.}

  \item{mode}{a character string specifying the type of
  penalty parameter.  If \code{"lambda"}, \code{lambda}
  gives penalty parameter directly.  If \code{"fraction"},
  the smallest value of the penalty parameter that sets all
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

  \item{use.correction}{currently ignored.  Small sample
  correction factors may be added in the future.}

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

  \item{ncores}{a positive integer giving the number of
  processor cores to be used for parallel computing (the
  default is 1 for no parallelization).  If this is set to
  \code{NA}, all available processor cores are used.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{model}{a logical indicating whether the data
  \code{x} and \code{y} should be added to the return
  object.  If \code{intercept} is \code{TRUE}, a column of
  ones is added to \code{x} to account for the intercept.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"sparseLTS"} with the following
  components:

  \item{best}{an integer vector giving the best subset of
  \eqn{h} observations found and used for computing the raw
  estimates.}

  \item{objective}{numeric; the value of the sparse LTS
  objective function, i.e., the \eqn{L_{1}}{L1} penalized
  sum of the \eqn{h} smallest squared residuals from the
  raw fit.}

  \item{coefficients}{a numeric vector of coefficient
  estimates of the reweighted fit (including the intercept
  if \code{intercept} is \code{TRUE}).}

  \item{fitted.values}{a numeric vector containing the
  fitted values of the response of the reweighted fit.}

  \item{residuals}{a numeric vector containing the
  residuals of the reweighted fit.}

  \item{center}{a numeric value giving the robust center
  estimate of the reweighted residuals.}

  \item{scale}{a numeric value giving the robust scale
  estimate of the reweighted residuals.}

  \item{lambda}{a numeric value giving the penalty
  parameter.}

  \item{intercept}{a logical indicating whether the model
  includes a constant term.}

  \item{alpha}{a numeric value giving the percentage of the
  residuals for which the \eqn{L_{1}}{L1} penalized sum of
  squares was minimized.}

  \item{quan}{the number \eqn{h} of observations used to
  compute the raw estimates.}

  \item{cnp2}{a numeric value giving the consistency factor
  applied to the scale estimate of the reweighted
  residuals.}

  \item{weights}{an integer vector containing binary
  weights that indicate outliers, i.e., the weights are
  \eqn{1} for observations with reasonably small reweighted
  residuals and \eqn{0} for observations with large
  reweighted residuals.}

  \item{df}{an integer giving the degrees of freedom of the
  obtained reweighted model fit, i.e., the number of
  nonzero coefficient estimates.}

  \item{raw.coefficients}{a numeric vector of coefficient
  estimates of the raw fit (including the intercept if
  \code{intercept} is \code{TRUE}).}

  \item{raw.residuals}{a numeric vector containing the
  residuals of the raw fit.}

  \item{raw.center}{a numeric value giving the robust
  center estimate of the raw residuals.}

  \item{raw.scale}{a numeric value giving the robust scale
  estimate of the raw residuals.}

  \item{raw.cnp2}{a numeric value giving the consistency
  factor applied to the scale estimate of the raw
  residuals.}

  \item{raw.weights}{an integer vector containing binary
  weights that indicate outliers of the raw fit, i.e., the
  weights used for the reweighted fit.}

  \item{x}{the predictor matrix (if \code{model} is
  \code{TRUE}).}

  \item{y}{the response variable (if \code{model} is
  \code{TRUE}).}

  \item{call}{the matched function call.}
}
\description{
  Compute least trimmed squares regression with an
  \eqn{L_{1}}{L1} penalty on the regression coefficients,
  which allows for sparse model estimates.
}
\note{
  Package \pkg{robustHD} has a built-in back end for sparse
  least trimmed squares using the C++ library Armadillo.
  Another back end is available through package
  \pkg{sparseLTSEigen}, which uses the C++ library Eigen.
  The latter is faster, but not available on all platforms.
  For instance, \pkg{sparseLTSEigen} currently does not
  work on 32-bit \R for Windows.  In addition, there is
  currently no binary package for OS X available on CRAN
  due to problems with the PowerPC architecture.
  Nevertheless, OS X users with Intel machines can install
  \pkg{RcppEigen} and \pkg{sparseLTSEigen} from source if
  the standard \R developer tools are installed.

  For both C++ back ends, parallel computing is implemented
  via OpenMP (\url{http://openmp.org/}).
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

## fit sparse LTS model
fit <- sparseLTS(x, y, lambda = 0.05, mode = "fraction")
coef(fit, zeros = FALSE)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{sparseLTSGrid}},
  \code{\link{coef.sparseLTS}},
  \code{\link{fitted.sparseLTS}},
  \code{\link{plot.sparseLTS}},
  \code{\link{predict.sparseLTS}},
  \code{\link{residuals.sparseLTS}},
  \code{\link{weights.sparseLTS}},
  \code{\link[robustbase]{ltsReg}}
}
\keyword{regression}
\keyword{robust}

