\name{plot.seqModel}
\alias{plot.grplars}
\alias{plot.rlars}
\alias{plot.seqModel}
\alias{plot.sparseLTSGrid}
\alias{plot.tslars}
\alias{plot.tslarsP}
\title{Plot a sequence of regression models}
\usage{
  \method{plot}{seqModel} (x,
    method = c("coefficients", "crit"), ...)

  \method{plot}{tslars} (x, p,
    method = c("coefficients", "crit"), ...)
}
\arguments{
  \item{x}{the model fit to be plotted.}

  \item{p}{an integer giving the lag length for which to
  produce the plot (the default is to use the optimal lag
  length).}

  \item{method}{a character string specifying the type of
  plot.  Possible values are \code{"coefficients"} to plot
  the coefficients from the submodels via
  \code{\link{coefPlot}}, or \code{"crit"} to plot the
  values of the optimality criterion for the submodels via
  \code{\link{critPlot}}.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"trellis"} (see
  \code{\link[lattice]{xyplot}}).
}
\description{
  Produce a plot of the coefficients or values of the
  optimality criterion for a sequence of regression models,
  such as submodels along a robust or groupwise least angle
  regression sequence, or sparse least trimmed squares
  regression models for a grid of values for the penalty
  parameter.
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

## create plots
plot(fit, zeros = FALSE)
plot(fit, method = "crit")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{coefPlot}}, \code{\link{critPlot}},
  \code{\link{rlars}}, \code{\link{grplars}},
  \code{\link{rgrplars}}, \code{\link{tslarsP}},
  \code{\link{rtslarsP}}, \code{\link{tslars}},
  \code{\link{rtslars}}, \code{\link{sparseLTSGrid}}
}
\keyword{hplot}

