\name{AIC.sparseLTS}
\alias{AIC.sparseLTS}
\alias{AIC.sparseLTSGrid}
\alias{BIC.sparseLTS}
\alias{BIC.sparseLTSGrid}
\title{Information criteria for sparse LTS regression models}
\usage{
  \method{AIC}{sparseLTS} (object, ...,
    fit = c("reweighted", "raw", "both"), k = 2)

  \method{BIC}{sparseLTS} (object, ...)

  \method{AIC}{sparseLTSGrid} (object, ...,
    fit = c("reweighted", "raw", "both"), k = 2)

  \method{BIC}{sparseLTSGrid} (object, ...)
}
\arguments{
  \item{object}{the model fit for which to compute the
  information criterion.}

  \item{\dots}{for the \code{BIC} method, additional
  arguments to be passed down to the \code{AIC} method.
  For the \code{AIC} method, additional arguments are
  currently ignored.}

  \item{fit}{a character string specifying for which fit to
  compute the information criterion.  Possible values are
  \code{"reweighted"} (the default) for the information
  criterion of the reweighted fit, \code{"raw"} for the
  information criterion of the raw fit, or \code{"both"}
  for the information criteria of both fits.}

  \item{k}{a numeric value giving the penalty per parameter
  to be used.  The default is to use \eqn{2} as in the
  classical definition of the AIC.}
}
\value{
  A numeric vector giving the information criteria for the
  requested fits.
}
\description{
  Compute the Akaike or Bayes information criterion for
  sparse least trimmed squares regression models based on
  the robust residual scale estimate.
}
\details{
  The information criteria are computed as \eqn{n (\log(2
  \pi) + 1 + \log(\hat{\sigma}^2)) + df k}{n (log(2 pi) + 1
  + log(sigma^2)) + df k}, where \eqn{n} denotes the number
  of observations, \eqn{\hat{\sigma}}{sigma} is the robust
  residual scale estimate, \eqn{df} is the number of
  nonzero coefficient estimates, and \eqn{k} is penalty per
  parameter.  The usual definition of the AIC uses \eqn{k =
  2}, whereas the BIC uses \eqn{k = \log(n)}{k = log(n)}.
  Consequently, the former is used as the default penalty
  of the \code{AIC} method, whereas the \code{BIC} method
  calls the \code{AIC} method with the latter penalty.
}
\note{
  Computing information criteria for several objects
  supplied via the \code{\dots} argument (as for the
  default methods of \code{\link[stats]{AIC}} and
  \code{BIC}) is currently not implemented.
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

## fit sparse LTS model over a grid of values for lambda
frac <- seq(0.25, 0.05, by = -0.05)
fitGrid <- sparseLTSGrid(x, y, lambda = frac, mode = "fraction")

## compute AIC and BIC
AIC(fitGrid)
BIC(fitGrid)
}
\author{
  Andreas Alfons
}
\references{
  Akaike, H. (1970) Statistical predictor identification.
  \emph{Annals of the Institute of Statistical
  Mathematics}, \bold{22}(2), 203--217.

  Schwarz, G. (1978) Estimating the dimension of a model.
  \emph{The Annals of Statistics}, \bold{6}(2), 461--464.
}
\seealso{
  \code{\link[stats]{AIC}}, \code{\link{sparseLTS}},
  \code{\link{sparseLTSGrid}}
}
\keyword{regression}

