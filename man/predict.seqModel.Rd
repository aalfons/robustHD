\name{predict.seqModel}
\alias{predict.grplars}
\alias{predict.optSeqModel}
\alias{predict.rlars}
\alias{predict.seqModel}
\alias{predict.tslars}
\alias{predict.tslarsP}
\title{Predict from a sequence of regression models}
\usage{
  \method{predict}{seqModel} (object, newdata, s = NA, ...)

  \method{predict}{optSeqModel} (object, newdata, ...)

  \method{predict}{tslarsP} (object, newdata, ...)

  \method{predict}{tslars} (object, newdata, p, ...)
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

  \item{p}{an integer giving the lag length for which to
  make predictions (the default is to use the optimal lag
  length).}

  \item{s}{an integer vector giving the steps of the
  submodels for which to make predictions (the default is
  to use the optimal submodel).}

  \item{\dots}{for the \code{"tslars"} method, additional
  arguments to be passed down to the \code{"tslarsP"}
  method.  For the other methods, additional arguments to
  be passed down to the respective method of
  \code{\link[=coef.seqModel]{coef}}.}
}
\value{
  If only one submodel is requested, a numeric vector
  containing the corresponding predicted values.

  If multiple submodels are requested, a numeric matrix in
  which each column contains the predicted values from the
  corresponding submodel.
}
\description{
  Make predictions from a sequence of regression models
  such as submodels along a robust or groupwise least angle
  regression sequence.  For autoregressive time series
  models with exogenous inputs, \eqn{h}-step ahead
  forecasts are performed.
}
\details{
  For \code{predict.seqModel}, the \code{newdata} argument
  defaults to the matrix of predictors used to fit the
  model such that the fitted values are computed.

  For autoregressive time series models with exogenous
  inputs with forecast horizon \eqn{h}, the \eqn{h} most
  recent observations of the predictors are omitted from
  fitting the model since there are no corresponding values
  for the response.  Hence the \code{newdata} argument for
  \code{predict.tslarsP} and \code{predict.tslars} defaults
  to those \eqn{h} observations of the predictors.
}
\examples{
## generate data
# example is not high-dimensional to keep computation time low
set.seed(1234)  # for reproducibility
n <- 100  # number of observations
p <- 25   # number of variables
beta <- rep.int(c(1, 0), c(5, p-5))  # coefficients
sigma <- 0.5      # controls signal-to-noise ratio
epsilon <- 0.1    # contamination level
x <- replicate(p, rnorm(n))     # predictor matrix
e <- rnorm(n)                   # error terms
i <- 1:ceiling(epsilon*n)       # observations to be contaminated
e[i] <- e[i] + 5                # vertical outliers
y <- c(x \%*\% beta + sigma * e)  # response
x[i,] <- x[i,] + 5              # bad leverage points

## fit robust LARS model
fit <- rlars(x, y)

## compute fitted values via predict method
predict(fit)
predict(fit, s = 1:5)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[stats]{predict}}, \code{\link{rlars}},
  \code{\link{grplars}}, \code{\link{rgrplars}},
  \code{\link{tslarsP}}, \code{\link{rtslarsP}},
  \code{\link{tslars}}, \code{\link{rtslars}}
}
\keyword{regression}

