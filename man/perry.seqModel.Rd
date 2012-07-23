\name{perry.seqModel}
\alias{perry.rlars}
\alias{perry.seqModel}
\title{Resampling-based prediction error for a sequential regression model}
\usage{
  \method{perry}{seqModel} (object, splits = foldControl(),
    cost, seed = NULL, ...)
}
\arguments{
  \item{object}{the model fit for which to estimate the
  prediction error.}

  \item{splits}{an object of class \code{"cvFolds"} (as
  returned by \code{\link[perry]{cvFolds}}) or a control
  object of class \code{"foldControl"} (see
  \code{\link[perry]{foldControl}}) defining the folds of
  the data for (repeated) \eqn{K}-fold cross-validation, an
  object of class \code{"randomSplits"} (as returned by
  \code{\link[perry]{randomSplits}}) or a control object of
  class \code{"splitControl"} (see
  \code{\link[perry]{splitControl}}) defining random data
  splits, or an object of class \code{"bootSamples"} (as
  returned by \code{\link[perry]{bootSamples}}) or a
  control object of class \code{"bootControl"} (see
  \code{\link[perry]{bootControl}}) defining bootstrap
  samples.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect vectors to be passed as its first two
  arguments, the first corresponding to the observed values
  of the response and the second to the predicted values,
  and must return a non-negative scalar value.  The default
  is to use the root mean squared prediction error for
  non-robust models and the root trimmed mean squared
  prediction error for robust models (see
  \code{\link[perry]{cost}}).}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed to the
  prediction loss function \code{cost}.}
}
\value{
  An object of class \code{"perry"} with the following
  components:

  \item{pe}{a numeric value giving the estimated prediction
  error.  In case of more than one replication, this gives
  the average value over all replications.}

  \item{se}{a numeric value giving the estimated standard
  error of the prediction loss.}

  \item{reps}{a numeric matrix with one column that
  contains the estimated prediction errors from all
  replications.  This is only returned in case of more than
  one replication.}

  \item{splits}{an object giving the data splits used to
  estimate the prediction error.}

  \item{y}{the response.}

  \item{yHat}{a list containing the predicted values from
  all replications.}

  \item{seed}{the seed of the random number generator
  before estimation of the prediction error.}

  \item{call}{the matched function call.}
}
\description{
  Estimate the prediction error of a previously fit
  sequential regression model such as a robust least angle
  regression model via (repeated) \eqn{K}-fold
  cross-validation, (repeated) random splitting (also known
  as random subsampling or Monte Carlo cross-validation),
  or the bootstrap.  In each iteration, the optimal model
  is thereby selected from the training data and used to
  make predictions for the test data.
}
\note{
  Users may prefer the wrapper functions
  \code{\link[perry]{repCV}}, \code{\link[perry]{repRS}}
  and \code{\link[perry]{bootPE}}.
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
repCV(fit)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{rlars}}, \code{\link{predict.seqModel}},
  \code{\link[perry]{perry}}, \code{\link[perry]{cost}}
}
\keyword{utilities}

