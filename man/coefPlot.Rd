\name{coefPlot}
\alias{coefPlot}
\alias{coefPlot.rlars}
\alias{coefPlot.seqModel}
\alias{coefPlot.sparseLTSGrid}
\title{Coefficient plot of a sequence of regression models}
\usage{
  coefPlot(x, ...)

  \method{coefPlot}{seqModel} (x,
    abscissa = c("step", "df"), zeros = FALSE, grid = TRUE,
    labels, pos = 4, offset = 0.5, ...)

  \method{coefPlot}{sparseLTSGrid} (x,
    fit = c("reweighted", "raw", "both"),
    abscissa = c("step", "df"), zeros = FALSE, grid = TRUE,
    labels, pos = 4, offset = 0.5, ...)
}
\arguments{
  \item{x}{the model fit to be plotted.}

  \item{fit}{a character string specifying for which
  estimator to produce the plot.  Possible values are
  \code{"reweighted"} (the default) for the reweighted
  fits, \code{"raw"} for the raw fits, or \code{"both"} for
  both estimators.}

  \item{abscissa}{a character string specifying what to
  plot on the \eqn{x}-axis.  Possible values are
  \code{"step"} for the step number (the default), or
  \code{"df"} for the degrees of freedom.}

  \item{zeros}{a logical indicating whether predictors that
  never enter the model and thus have zero coefficients
  should be included in the plot (\code{TRUE}) or omitted
  (\code{FALSE}, the default).  This is useful if the
  number of predictors is much larger than the number of
  observations, in which case many coefficients are never
  nonzero.}

  \item{grid}{a logical indicating whether vertical grid
  lines should be drawn at each step.}

  \item{labels}{an optional character vector containing
  labels for the predictors.}

  \item{pos}{an integer position specifier for the labels.
  Possible values are 1, 2, 3 and 4, respectively
  indicating positions below, to the left of, above and to
  the right of the corresponding coefficient values from
  the last step.}

  \item{offset}{a numeric value giving the offset of the
  labels from the corresponding coefficient values from the
  last step (in fractions of a character width).}

  \item{\dots}{for the generic function, additional
  arguments to be passed down to methods.  For the
  \code{"seqModel"} and \code{"sparseLTSGrid"} methods,
  additional arguments to be passed down to
  \code{\link[lattice]{xyplot}}.}
}
\value{
  An object of class \code{"trellis"} (see
  \code{\link[lattice]{xyplot}}).
}
\description{
  Produce a plot of the coefficients from a sequence of
  regression models, such as submodels along a robust least
  angle regression sequence, or sparse least trimmed
  squares regression models for a grid of values for the
  penalty parameter.
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
coefPlot(fit, zeros = FALSE)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[lattice]{xyplot}}, \code{\link{rlars}},
  \code{\link{sparseLTSGrid}}
}
\keyword{hplot}

