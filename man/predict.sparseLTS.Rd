\name{predict.sparseLTS}
\alias{predict.sparseLTS}
\title{Predict from sparse LTS regression models}
\usage{
  \method{predict}{sparseLTS} (object, newdata, s = NA,
    fit = c("reweighted", "raw", "both"), ...)
}
\arguments{
  \item{object}{the model fit from which to make
  predictions.}

  \item{newdata}{new data for the predictors.  If the model
  fit was computed with the formula method, this should be
  a data frame from which to extract the predictor
  variables.  Otherwise this should be a matrix containing
  the same variables as the predictor matrix used to fit
  the model (possibly including a column of ones in case of
  a model with intercept).}

  \item{s}{an integer vector giving the indices of the
  models for which to make predictions.  If \code{fit} is
  \code{"both"}, this can be a list with two components,
  with the first component giving the indices of the
  reweighted fits and the second the indices of the raw
  fits.  The default is to use the optimal model for each
  of the requested estimators.  Note that the optimal
  models may not correspond to the same value of the
  penalty parameter for the reweighted and the raw
  estimator.}

  \item{fit}{a character string specifying for which fit to
  make predictions.  Possible values are
  \code{"reweighted"} (the default) for predicting values
  from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit, or \code{"both"} for predicting
  values from both fits.}

  \item{\dots}{currently ignored.}
}
\value{
  If predictions for only one model are requested, they are
  returned in the form of a numeric vector.

  Otherwise a numeric matrix is returned in which each
  column contains the predicted values from the
  corresponding model.
}
\description{
  Make predictions from sparse least trimmed squares
  regression models.
}
\details{
  The \code{newdata} argument defaults to the matrix of
  predictors used to fit the model such that the fitted
  values are computed.
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
# compute fitted values via predict method
predict(fit)
predict(fit, fit = "both")

## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.25, 0.05, by = -0.05)
fitGrid <- sparseLTS(x, y, lambda = frac, mode = "fraction")
# compute fitted values via predict method
predict(fitGrid)
predict(fitGrid, fit = "both")
predict(fitGrid, s = NULL)
predict(fitGrid, fit = "both", s = NULL)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[stats]{predict}}, \code{\link{sparseLTS}}
}
\keyword{regression}

