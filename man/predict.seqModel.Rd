\name{predict.seqModel}
\alias{predict.rlars}
\alias{predict.seqModel}
\alias{predict.sparseLTS}
\title{Predict from a sequence of regression models}
\usage{
  \method{predict}{seqModel} (object, newdata, s = NA, ...)

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
  the model (including a column of ones to account for the
  intercept).}

  \item{s}{for the \code{"seqModel"} method, an integer
  vector giving the steps of the submodels for which to
  make predictions (the default is to use the optimal
  submodel).  For the \code{"sparseLTS"} method, an integer
  vector giving the indices of the models for which to make
  predictions.  If \code{fit} is \code{"both"}, this can be
  a list with two components, with the first component
  giving the indices of the reweighted fits and the second
  the indices of the raw fits.  The default is to use the
  optimal model for each of the requested estimators.  Note
  that the optimal models may not correspond to the same
  value of the penalty parameter for the reweighted and the
  raw estimator.}

  \item{fit}{a character string specifying for which fit to
  make predictions.  Possible values are
  \code{"reweighted"} (the default) for predicting values
  from the reweighted fit, \code{"raw"} for predicting
  values from the raw fit, or \code{"both"} for predicting
  values from both fits.}

  \item{\dots}{additional arguments to be passed down to
  the respective method of
  \code{\link[=coef.seqModel]{coef}}.}
}
\value{
  A numeric vector or matrix containing the requested
  predicted values.
}
\description{
  Make predictions from a sequence of regression models,
  such as submodels along a robust least angle regression
  sequence, or sparse least trimmed squares regression
  models for a grid of values for the penalty parameter.
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


## robust LARS
# fit model
fitRlars <- rlars(x, y, sMax = 25)
# extract coefficients
predict(fitRlars)
predict(fitRlars, s = 1:5)


## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.25, 0.05, by = -0.05)
fitSparseLTS <- sparseLTS(x, y, lambda = frac, mode = "fraction")
# extract coefficients
predict(fitSparseLTS)
predict(fitSparseLTS, fit = "both")
predict(fitSparseLTS, s = NULL)
predict(fitSparseLTS, fit = "both", s = NULL)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[stats]{predict}}, \code{\link{rlars}},
  \code{\link{sparseLTS}}
}
\keyword{regression}

