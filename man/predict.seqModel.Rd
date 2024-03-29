% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.seqModel}
\alias{predict.seqModel}
\alias{predict.rlars}
\alias{predict.grplars}
\alias{predict.tslarsP}
\alias{predict.tslars}
\alias{predict.sparseLTS}
\title{Predict from a sequence of regression models}
\usage{
\method{predict}{seqModel}(object, newdata, s = NA, ...)

\method{predict}{tslarsP}(object, newdata, ...)

\method{predict}{tslars}(object, newdata, p, ...)

\method{predict}{sparseLTS}(object, newdata, s = NA, fit = c("reweighted", "raw", "both"), ...)
}
\arguments{
\item{object}{the model fit from which to make predictions.}

\item{newdata}{new data for the predictors.  If the model fit was computed
with the formula method, this should be a data frame from which to extract
the predictor variables.  Otherwise this should be a matrix containing the
same variables as the predictor matrix used to fit the model (including a
column of ones to account for the intercept).}

\item{s}{for the \code{"seqModel"} method, an integer vector giving the
steps of the submodels for which to make predictions (the default is to use
the optimal submodel).  For the \code{"sparseLTS"} method, an integer vector
giving the indices of the models for which to make predictions.  If
\code{fit} is \code{"both"}, this can be a list with two components, with
the first component giving the indices of the reweighted fits and the second
the indices of the raw fits.  The default is to use the optimal model for
each of the requested estimators.  Note that the optimal models may not
correspond to the same value of the penalty parameter for the reweighted and
the raw estimator.}

\item{\dots}{for the \code{"tslars"} method, additional arguments to be
passed down to the \code{"tslarsP"} method.  For the other methods,
additional arguments to be passed down to the respective method of
\code{\link[=coef.seqModel]{coef}}.}

\item{p}{an integer giving the lag length for which to make predictions
(the default is to use the optimal lag length).}

\item{fit}{a character string specifying for which fit to make
predictions.  Possible values are \code{"reweighted"} (the default) for
predicting values from the reweighted fit, \code{"raw"} for predicting
values from the raw fit, or \code{"both"} for predicting values from both
fits.}
}
\value{
A numeric vector or matrix containing the requested predicted values.
}
\description{
Make predictions from a sequence of regression models, such as submodels
along a robust or groupwise least angle regression sequence, or sparse least
trimmed squares regression models for a grid of values for the penalty
parameter.  For autoregressive time series models with exogenous inputs,
\eqn{h}-step ahead forecasts are performed.
}
\details{
The \code{newdata} argument defaults to the matrix of predictors used to fit
the model such that the fitted values are computed.

For autoregressive time series models with exogenous inputs with forecast
horizon \eqn{h}, the \eqn{h} most recent observations of the predictors are
omitted from fitting the model since there are no corresponding values for
the response.  Hence the \code{newdata} argument for \code{predict.tslarsP}
and \code{predict.tslars} defaults to those \eqn{h} observations of the
predictors.
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
fitRlars <- rlars(x, y, sMax = 10)
# compute fitted values via predict method
predict(fitRlars)
head(predict(fitRlars, s = 1:5))


## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.2, 0.05, by = -0.05)
fitSparseLTS <- sparseLTS(x, y, lambda = frac, mode = "fraction")
# compute fitted values via predict method
predict(fitSparseLTS)
head(predict(fitSparseLTS, fit = "both"))
head(predict(fitSparseLTS, s = NULL))
head(predict(fitSparseLTS, fit = "both", s = NULL))
}
\seealso{
\code{\link[stats]{predict}}, \code{\link{rlars}},
\code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
\code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
\code{\link{sparseLTS}}
}
\author{
Andreas Alfons
}
\keyword{regression}
