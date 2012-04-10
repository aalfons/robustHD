\name{residuals.seqModel}
\alias{residuals.rlars}
\alias{residuals.seqModel}
\title{Extract residuals from a sequence of regression models}
\usage{
  \method{residuals}{seqModel} (object, s, ...)
}
\arguments{
  \item{object}{the model fit from which to extract
  residuals.}

  \item{s}{an integer vector giving the steps of the
  submodels for which to extract the residuals (the default
  is to use the optimal submodel).}

  \item{\dots}{additional arguments are currently ignored.}
}
\value{
  If only one submodel is requested, a numeric vector
  containing the corresponding residuals.

  If multiple submodels are requested, a numeric matrix in
  which each column contains the residuals of the
  corresponding submodel.
}
\description{
  Extract residuals from a sequence of regression models
  such as submodels along a robust least angle regression
  sequence.
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
  \code{\link[stats]{residuals}}, \code{\link{rlars}}
}
\keyword{regression}

