\name{cvRlars}
\alias{cvRlars}
\alias{cvRlars.default}
\alias{cvRlars.formula}
\title{Cross-validation along a robust least angle regression sequence}
\usage{
  cvRlars(x, ...)

  \method{cvRlars}{formula} (formula, data, ...)

  \method{cvRlars}{default} (x, y, cost = rtmspe, K = 5,
    R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, selectBest = c("min", "hastie"),
    seFactor = 1, active, s = NULL, regFun = lmrob,
    regArgs = list(), seed = NULL, ...)
}
\arguments{
  \item{formula}{a formula describing the full model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{cvRlars} is called.}

  \item{x}{a matrix or data frame containing the candidate
  predictors.}

  \item{y}{a numeric vector containing the response.}

  \item{cost}{a robust cost function measuring prediction
  loss.  It should expect vectors to be passed as its first
  two arguments, the first corresponding to the observed
  values of the response and the second to the predicted
  values, and must return a non-negative scalar value.  The
  default is to use the root trimmed mean squared
  prediction error (see \code{\link[cvTools]{cost}}).}

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

  \item{selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter selects the most parsimonious model
  whose prediction error is no larger than \code{seFactor}
  standard errors above the prediction error of the best
  overall model.}

  \item{seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{selectBest} is
  \code{"min"}.}

  \item{active}{an integer vector containing the sequence
  of predictor groups (as returned by
  \code{\link{rlars}}).}

  \item{s}{an integer vector giving the steps of the
  submodels for which to estimate the prediction errors
  (the default is to use all steps along the sequence as
  long as there are twice as many observations as
  predictors).}

  \item{regFun}{a function to compute robust linear
  regressions (defaults to
  \code{\link[robustbase]{lmrob}}).}

  \item{regArgs}{a list of arguments to be passed to
  \code{regFun}.}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed to the
  prediction loss function \code{cost}.}
}
\value{
  An object of class \code{"cvSeqModel"} (which inherits
  from class \code{"cvSelect"}) with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds used in
  cross-validation.}

  \item{R}{an integer giving the number of replications
  used in cross-validation.}

  \item{best}{an integer giving the index of the submodel
  with the best prediction performance.}

  \item{cv}{a data frame containing the estimated
  prediction errors for the submodels.  For repeated
  cross-validation, those are average values over all
  replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for the submodels.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{reps}{a data frame containing the estimated
  prediction errors for the submodels from all
  replications.  This is only returned for repeated
  cross-validation.}

  \item{call}{the matched function call.}
}
\description{
  Estimate the prediction error of submodels along a robust
  least angle regression sequence via (repeated)
  \eqn{K}-fold cross-validation.
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

## obtain robust LARS sequence
active <- rlars(x, y, fit = FALSE)

## evaluate models along sequence
cv <- cvRlars(x, y, active = active, selectBest = "hastie", 
    includeSE = TRUE)
cv
dotplot(cv)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{repCV.rlars}}, \code{\link{rlars}},
  \code{\link{predict.rlars}},
  \code{\link[cvTools]{cvFolds}},
  \code{\link[cvTools]{cost}}
}
\keyword{robust}
\keyword{utilities}

