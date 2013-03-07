\name{residuals.seqModel}
\alias{residuals.rlars}
\alias{residuals.seqModel}
\alias{residuals.sparseLTS}
\title{Extract residuals from a sequence of regression models}
\usage{
  \method{residuals}{seqModel} (object, s = NA,
    standardized = FALSE, drop = !is.null(s), ...)

  \method{residuals}{sparseLTS} (object, s = NA,
    fit = c("reweighted", "raw", "both"),
    standardized = FALSE, drop = !is.null(s), ...)
}
\arguments{
  \item{object}{the model fit from which to extract
  residuals.}

  \item{s}{for the \code{"seqModel"} method, an integer
  vector giving the steps of the submodels for which to
  extract the residuals (the default is to use the optimal
  submodel).  For the \code{"sparseLTS"} method, an integer
  vector giving the indices of the models for which to
  extract residuals.  If \code{fit} is \code{"both"}, this
  can be a list with two components, with the first
  component giving the indices of the reweighted fits and
  the second the indices of the raw fits.  The default is
  to use the optimal model for each of the requested
  estimators.  Note that the optimal models may not
  correspond to the same value of the penalty parameter for
  the reweighted and the raw estimator.}

  \item{fit}{a character string specifying which residuals
  to extract. Possible values are \code{"reweighted"} (the
  default) for the residuals from the reweighted estimator,
  \code{"raw"} for the residuals from the raw estimator, or
  \code{"both"} for the residuals from both estimators.}

  \item{standardized}{a logical indicating whether the
  residuals should be standardized (the default is
  \code{FALSE}).}

  \item{drop}{a logical indicating whether to reduce the
  dimension to a vector in case of only one step.}

  \item{\dots}{additional arguments are currently ignored.}
}
\value{
  A numeric vector or matrix containing the requested
  residuals.
}
\description{
  Extract residuals from a sequence of regression models,
  such as submodels along a robust least angle regression
  sequence, or sparse least trimmed squares regression
  models for a grid of values for the penalty parameter.
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


## robust LARS
# fit model
fitRlars <- rlars(x, y, sMax = 25)
# extract coefficients
residuals(fitRlars)
residuals(fitRlars, s = 0:4)


## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.25, 0.05, by = -0.05)
fitSparseLTS <- sparseLTS(x, y, lambda = frac, mode = "fraction")
# extract coefficients
residuals(fitSparseLTS)
residuals(fitSparseLTS, fit = "both")
residuals(fitSparseLTS, s = NULL)
residuals(fitSparseLTS, fit = "both", s = NULL)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[stats]{residuals}}, \code{\link{rlars}},
  \code{\link{sparseLTS}}
}
\keyword{regression}

