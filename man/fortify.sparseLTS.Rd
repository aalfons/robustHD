\name{fortify.sparseLTS}
\alias{fortify.sparseLTS}
\alias{fortify.sparseLTSGrid}
\title{Convert a sparse LTS regression model into a data frame for plotting}
\usage{
  \method{fortify}{sparseLTS} (model, data,
    fit = c("reweighted", "raw", "both"), ...)

  \method{fortify}{sparseLTSGrid} (model, data, s,
    fit = c("reweighted", "raw", "both"), ...)
}
\arguments{
  \item{model}{the model fit to be converted.}

  \item{data}{currently ignored.}

  \item{s}{an integer vector giving the indices of the
  models to be converted.  The default is to use the
  optimal model for each of the requested fits.}

  \item{fit}{a character string specifying which fit to
  convert.  Possible values are \code{"reweighted"} (the
  default) to convert the reweighted fit, \code{"raw"} to
  convert the raw fit, or \code{"both"} to convert both
  fits.}

  \item{\dots}{currently ignored.}
}
\value{
  A data frame containing the columns listed below, as well
  as additional information stored in the attributes
  \code{"qqLine"} (intercepts and slopes of the respective
  reference lines to be displayed in residual Q-Q plots),
  \code{"q"} (quantiles of the Mahalanobis distribution
  used as cutoff points for detecting leverage points), and
  \code{"facets"} (default faceting formula for the
  diagnostic plots).

  \item{step}{the indices of the models (only returned if
  more than one model is requested).}

  \item{fit}{the model fits (only returned if both the
  reweighted and raw fit are requested).}

  \item{index}{the indices of the observations.}

  \item{fitted}{the fitted values.}

  \item{residual}{the standardized residuals.}

  \item{theoretical}{the corresponding theoretical
  quantiles from the standard normal distribution.}

  \item{qqd}{the absolute distances from a reference line
  through the first and third sample and theoretical
  quartiles.}

  \item{rd}{the robust Mahalanobis distances computed via
  the MCD (see \code{\link[robustbase]{covMcd}}).}

  \item{xyd}{the pairwise maxima of the absolute values of
  the standardized residuals and the robust Mahalanobis
  distances, divided by the respective other outlier
  detection cutoff point.}

  \item{weight}{the weights indicating regression
  outliers.}

  \item{leverage}{logicals indicating leverage points
  (i.e., outliers in the predictor space).}

  \item{classification}{a factor with levels
  \code{"outlier"} (regression outliers) and \code{"good"}
  (data points following the model).}
}
\description{
  Supplement the fitted values and residuals of a sparse
  least trimmed squared model with other useful information
  for diagnostic plots.
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

## sparse LTS
# fit model
fit <- sparseLTS(x, y, lambda = 0.05, mode = "fraction")
# create plot
fortify(fit)
fortify(fit, fit = "both")

## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.25, 0.05, by = -0.05)
fitGrid <- sparseLTSGrid(x, y, lambda = frac, mode = "fraction")
# create plot
fortify(fitGrid)
fortify(fitGrid, fit = "both")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[ggplot2]{fortify}},
  \code{\link{diagnosticPlot}}, \code{\link{sparseLTS}},
  \code{\link{sparseLTSGrid}}
}
\keyword{utilities}

