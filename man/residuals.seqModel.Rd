\name{residuals.seqModel}
\alias{residuals.grplars}
\alias{residuals.optSeqModel}
\alias{residuals.rlars}
\alias{residuals.seqModel}
\alias{residuals.tslars}
\alias{residuals.tslarsP}
\title{Extract residuals from a sequence of regression models}
\usage{
  \method{residuals}{seqModel} (object, s = NA,
    drop = !is.null(s), ...)

  \method{residuals}{optSeqModel} (object, ...)

  \method{residuals}{tslars} (object, p, ...)
}
\arguments{
  \item{object}{the model fit from which to extract
  residuals.}

  \item{p}{an integer giving the lag length for which to
  extract residuals (the default is to use the optimal lag
  length).}

  \item{s}{an integer vector giving the steps of the
  submodels for which to extract the residuals (the default
  is to use the optimal submodel).}

  \item{drop}{a logical indicating whether to reduce the
  dimension to a vector in case of only one step.}

  \item{\dots}{for the \code{"tslars"} method, additional
  arguments to be passed down to the \code{"seqModel"} or
  \code{"optSeqModel"} method.  For the other methods,
  additional arguments are currently ignored.}
}
\value{
  A numeric vector or matrix containing the requested
  residuals.
}
\description{
  Extract residuals from a sequence of regression models
  such as submodels along a robust or groupwise least angle
  regression sequence.
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

## extract residuals
residuals(fit)
residuals(fit, s = 1:5)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[stats]{residuals}}, \code{\link{rlars}},
  \code{\link{grplars}}, \code{\link{rgrplars}},
  \code{\link{tslarsP}}, \code{\link{rtslarsP}},
  \code{\link{tslars}}, \code{\link{rtslars}}
}
\keyword{regression}

