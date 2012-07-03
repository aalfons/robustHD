\name{critPlot}
\alias{critPlot}
\alias{critPlot.rlars}
\alias{critPlot.seqModel}
\alias{critPlot.sparseLTSGrid}
\title{Optimality criterion plot of a sequence of regression models}
\usage{
  critPlot(x, ...)

  \method{critPlot}{seqModel} (x, size = c(0.5, 2), ...)

  \method{critPlot}{sparseLTSGrid} (x,
    fit = c("reweighted", "raw", "both"), size = c(0.5, 2),
    ...)
}
\arguments{
  \item{x}{the model fit to be plotted.}

  \item{fit}{a character string specifying for which
  estimator to produce the plot.  Possible values are
  \code{"reweighted"} (the default) for the reweighted
  fits, \code{"raw"} for the raw fits, or \code{"both"} for
  both estimators.}

  \item{size}{a numeric vector of length two giving the
  line width and the point size, respectively.}

  \item{\dots}{for the generic function, additional
  arguments to be passed down to methods.  For the
  \code{"seqModel"} and \code{"sparseLTSGrid"} methods,
  additional arguments to be passed down to
  \code{\link[ggplot2]{geom_line}} and
  \code{\link[ggplot2]{geom_point}}.}
}
\value{
  An object of class \code{"ggplot"} (see
  \code{\link[ggplot2]{ggplot}}).
}
\description{
  Produce a plot of the values of the optimality criterion
  for a sequence of regression models, such as submodels
  along a robust least angle regression sequence, or sparse
  least trimmed squares regression models for a grid of
  values for the penalty parameter.
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

## create plot
critPlot(fit)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}},
  \code{\link{sparseLTSGrid}}
}
\keyword{hplot}

