\name{plot.seqModel}
\alias{plot.grplars}
\alias{plot.perrySeqModel}
\alias{plot.perrySparseLTS}
\alias{plot.rlars}
\alias{plot.seqModel}
\alias{plot.sparseLTS}
\alias{plot.tslars}
\alias{plot.tslarsP}
\title{Plot a sequence of regression models}
\usage{
  \method{plot}{seqModel} (x,
    method = c("coefficients", "crit", "diagnostic"), ...)

  \method{plot}{perrySeqModel} (x,
    method = c("crit", "diagnostic"), ...)

  \method{plot}{tslars} (x, p,
    method = c("coefficients", "crit", "diagnostic"), ...)

  \method{plot}{sparseLTS} (x,
    method = c("coefficients", "crit", "diagnostic"), ...)

  \method{plot}{perrySparseLTS} (x,
    method = c("crit", "diagnostic"), ...)
}
\arguments{
  \item{x}{the model fit to be plotted.}

  \item{p}{an integer giving the lag length for which to
  produce the plot (the default is to use the optimal lag
  length).}

  \item{method}{a character string specifying the type of
  plot.  Possible values are \code{"coefficients"} to plot
  the coefficients from the submodels via
  \code{\link{coefPlot}} (only for the \code{"seqModel"}
  and \code{"sparseLTS"} methods), \code{"crit"} to plot
  the values of the optimality criterion for the submodels
  via \code{\link{critPlot}}, or \code{"diagnostic"} for
  diagnostic plots via \code{\link{diagnosticPlot}}.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"ggplot"} (see
  \code{\link[ggplot2]{ggplot}}).
}
\description{
  Produce a plot of the coefficients, the values of the
  optimality criterion, or diagnostic plots for a sequence
  of regression models, such as submodels along a robust or
  groupwise least angle regression sequence, or sparse
  least trimmed squares regression models for a grid of
  values for the penalty parameter.
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
# create plots
plot(fitRlars, method = "coef")
plot(fitRlars, method = "crit")
plot(fitRlars, method = "diagnostic")


## sparse LTS over a grid of values for lambda
# fit model
frac <- seq(0.2, 0.05, by = -0.05)
fitSparseLTS <- sparseLTS(x, y, lambda = frac, mode = "fraction")
# create plots
plot(fitSparseLTS, method = "coef")
plot(fitSparseLTS, method = "crit")
plot(fitSparseLTS, method = "diagnostic")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{coefPlot}}, \code{\link{critPlot}},
  \code{\link{diagnosticPlot}}, \code{\link{rlars}},
  \code{\link{grplars}}, \code{\link{rgrplars}},
  \code{\link{tslarsP}}, \code{\link{rtslarsP}},
  \code{\link{tslars}}, \code{\link{rtslars}},
  \code{\link{sparseLTS}}
}
\keyword{hplot}

