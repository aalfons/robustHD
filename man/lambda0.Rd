\name{lambda0}
\alias{lambda0}
\title{Penalty parameter for sparse LTS regression}
\usage{
  lambda0(x, y, normalize = TRUE, intercept = TRUE,
    const = 2, prob = 0.95, tol = .Machine$double.eps^0.5,
    eps = .Machine$double.eps, ...)
}
\arguments{
  \item{x}{a numeric matrix containing the predictor
  variables.}

  \item{y}{a numeric vector containing the response
  variable.}

  \item{normalize}{a logical indicating whether the
  winsorized predictor variables should be normalized to
  have unit \eqn{L_{2}}{L2} norm (the default is
  \code{TRUE}).}

  \item{intercept}{a logical indicating whether a constant
  term should be included in the model (the default is
  \code{TRUE}).}

  \item{const}{numeric; tuning constant to be used in
  univariate winsorization (defaults to 2).}

  \item{prob}{numeric; probability for the quantile of the
  \eqn{\chi^{2}}{chi-squared} distribution to be used in
  bivariate winsorization (defaults to 0.95).}

  \item{tol}{a small positive numeric value used to
  determine singularity issues in the computation of
  correlation estimates for bivariate winsorization (see
  \code{\link{corHuber}}).}

  \item{eps}{a small positive numeric value used to
  determine whether the robust scale estimate of a variable
  is too small (an effective zero).}

  \item{\dots}{additional arguments to be passed to
  \code{\link[=standardize]{robStandardize}}.}
}
\value{
  A robust estimate of the smallest value of the penalty
  parameter for sparse LTS regression that sets all
  coefficients to zero.
}
\description{
  Use bivariate winsorization to estimate the smallest
  value of the penalty parameter for sparse least trimmed
  squares regression that sets all coefficients to zero.
}
\details{
  The estimation procedure is inspired by the calculation
  of the respective penalty parameter in the first step of
  the classical LARS algorithm. First, two-dimensional data
  blocks consisting of the response with each predictor
  variable are cleaned via bivariate winsorization.  For
  each block, the following computations are then
  performed.  If an intercept is included in the model, the
  cleaned response is centered and the corresponding
  cleaned predictor is centered and scaled to have unit
  norm.  Otherwise the variables are not centered, but the
  predictor is scaled to have unit norm.  Finally, the dot
  product of the response and the corresponding predictor
  is computed.  The largest absolute value of those dot
  products, rescaled to fit the parametrization of the
  sparse LTS definition, yields the estimate of the
  smallest penalty parameter that sets all coefficients to
  zero.
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

## estimate smallest value of the penalty parameter 
## that sets all coefficients to 0
lambda0(x, y)
}
\author{
  Andreas Alfons
}
\references{
  Alfons, A., Croux, C. and Gelper, S. (2013) Sparse least
  trimmed squares regression for analyzing high-dimensional
  large data sets. \emph{The Annals of Applied Statistics},
  \bold{7}(1), 226--248.

  Efron, B., Hastie, T., Johnstone, I. and Tibshirani, R.
  (2004) Least angle regression. \emph{The Annals of
  Statistics}, \bold{32}(2), 407--499.

  Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust
  linear model selection based on least angle regression.
  \emph{Journal of the American Statistical Association},
  \bold{102}(480), 1289--1299.
}
\seealso{
  \code{\link{sparseLTS}}, \code{\link{winsorize}}
}
\keyword{robust}
\keyword{utilities}

