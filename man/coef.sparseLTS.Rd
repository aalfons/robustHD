\name{coef.sparseLTS}
\alias{coef.sparseLTS}
\alias{coef.sparseLTSGrid}
\title{Extract coefficients from sparse LTS regression models}
\usage{
  \method{coef}{sparseLTS} (object,
    fit = c("reweighted", "raw", "both"), zeros = TRUE,
    ...)

  \method{coef}{sparseLTSGrid} (object, s,
    fit = c("reweighted", "raw", "both"), zeros = TRUE,
    ...)
}
\arguments{
  \item{object}{the model fit from which to extract
  coefficients.}

  \item{s}{an integer vector giving the indices of the
  models for which to extract coefficients.  If \code{fit}
  is \code{"both"}, this can be a list with two components,
  with the first component giving the indices of the
  reweighted fits and the second the indices of the raw
  fits.  The default is to use the optimal model for each
  of the requested estimators.  Note that the optimal
  models may not correspond to the same value of the
  penalty parameter for the reweighted and the raw
  estimator.}

  \item{fit}{a character string specifying which
  coefficients to extract. Possible values are
  \code{"reweighted"} (the default) for the coefficients
  from the reweighted estimator, \code{"raw"} for the
  coefficients from the raw estimator, or \code{"both"} for
  the coefficients from both estimators.}

  \item{zeros}{a logical indicating whether to keep zero
  coefficients (\code{TRUE}, the default) or to drop them
  (\code{FALSE}).}

  \item{\dots}{currently ignored.}
}
\value{
  If coefficients for only one model are requested, they
  are returned in the form of a numeric vector.

  Otherwise a numeric matrix is returned in which each
  column contains the coefficients of the corresponding
  model.
}
\description{
  Extract coefficients from sparse least trimmed squares
  regression models.
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

## sparse LTS
# fit model
fit <- sparseLTS(x, y, lambda = 0.05, mode = "fraction")
# extract coefficients
coef(fit, zeros = FALSE)
coef(fit, fit = "both", zeros = FALSE)

## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.25, 0.05, by = -0.05)
fitGrid <- sparseLTSGrid(x, y, lambda = frac, mode = "fraction")
# extract coefficients
coef(fitGrid, zeros = FALSE)
coef(fitGrid, fit = "both", zeros = FALSE)
coef(fitGrid, s = NULL, zeros = FALSE)
coef(fitGrid, fit = "both", s = NULL, zeros = FALSE)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[stats]{coef}}, \code{\link{sparseLTS}},
  \code{\link{sparseLTSGrid}}
}
\keyword{regression}

