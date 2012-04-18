\name{repCV.seqModel}
\alias{repCV.grplars}
\alias{repCV.rlars}
\alias{repCV.seqModel}
\title{Cross-validation for a sequential regression model}
\usage{
  \method{repCV}{seqModel} (object, cost, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, seed = NULL, ...)
}
\arguments{
  \item{object}{the model fit for which to estimate the
  prediction error.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect vectors to be passed as its first two
  arguments, the first corresponding to the observed values
  of the response and the second to the predicted values,
  and must return a non-negative scalar value.  The default
  is to use the root mean squared prediction error for
  non-robust models and the root trimmed mean squared
  prediction error for robust models (see
  \code{\link[cvTools]{cost}}).}

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

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed to the
  prediction loss function \code{cost}.}
}
\value{
  An object of class \code{"cv"} with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{cv}{a numeric value giving the estimated prediction
  error.  For repeated cross-validation, this gives the
  average value over all replications.}

  \item{se}{a numeric value giving the estimated standard
  error of the prediction loss.}

  \item{reps}{a numeric matrix with one column that
  contains the estimated prediction errors from all
  replications.  This is only returned for repeated
  cross-validation.}

  \item{seed}{the seed of the random number generator
  before cross-validation was performed.}

  \item{call}{the matched function call.}
}
\description{
  Perform (repeated) \eqn{K}-fold cross-validation to
  estimate the prediction error of a previously fit
  sequential regression model such as a robust or groupwise
  least angle regression model.  In each iteration of
  cross-validation, the optimal model is thereby selected
  from the training data and used to make predictions for
  the test data.
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

## fit and evaluate robust LARS model
fit <- rlars(x, y)
cv <- repCV(fit)
cv
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{rlars}}, \code{\link{grplars}},
  \code{\link{rgrplars}}, \code{\link{predict.seqModel}},
  \code{\link[cvTools]{cvFolds}},
  \code{\link[cvTools]{cost}}
}
\keyword{utilities}

