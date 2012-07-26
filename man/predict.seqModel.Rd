\name{predict.seqModel}
\alias{predict.optSeqModel}
\alias{predict.rlars}
\alias{predict.seqModel}
\title{Predict from a sequence of regression models}
\usage{
  \method{predict}{seqModel} (object, newdata, s = NA, ...)

  \method{predict}{optSeqModel} (object, newdata, ...)
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

  \item{s}{an integer vector giving the steps of the
  submodels for which to make predictions (the default is
  to use the optimal submodel).}

  \item{\dots}{additional arguments to be passed down to
  the respective method of
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
  such as submodels along a robust least angle regression
  sequence.
}
\details{
  For \code{predict.seqModel}, the \code{newdata} argument
  defaults to the matrix of predictors used to fit the
  model such that the fitted values are computed.
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
  \code{\link[stats]{predict}}, \code{\link{rlars}}
}
\keyword{regression}

