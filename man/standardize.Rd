\name{standardize}
\alias{robStandardize}
\alias{standardize}
\title{Data standardization}
\usage{
  standardize(x, centerFun = mean, scaleFun = sd)

  robStandardize(x, centerFun = median, scaleFun = mad,
    fallback = FALSE, eps = .Machine$double.eps, ...)
}
\arguments{
  \item{x}{a numeric vector, matrix or data frame to be
  standardized.}

  \item{centerFun}{a function to compute an estimate of the
  center of a variable (defaults to \code{\link{mean}}).}

  \item{scaleFun}{a function to compute an estimate of the
  scale of a variable (defaults to
  \code{\link[stats]{sd}}).}

  \item{fallback}{a logical indicating whether
  standardization with \code{\link{mean}} and
  \code{\link[stats]{sd}} should be performed as a fallback
  mode for variables whose robust scale estimate is too
  small.  This is useful, e.g., for data containing dummy
  variables.}

  \item{eps}{a small positive numeric value used to
  determine whether the robust scale estimate of a variable
  is too small (an effective zero).}

  \item{\dots}{currently ignored.}
}
\value{
  An object of the same type as the original data \code{x}
  containing the centered and scaled data.  The center and
  scale estimates of the original data are returned as
  attributes \code{"center"} and \code{"scale"},
  respectively.
}
\description{
  Standardize data with given functions for computing
  center and scale.
}
\details{
  \code{robStandardize} is a wrapper function for robust
  standardization, hence the default is to use
  \code{\link[stats]{median}} and \code{\link[stats]{mad}}.
}
\note{
  The implementation contains special cases for the
  typically used combinations
  \code{\link{mean}}/\code{\link[stats]{sd}} and
  \code{\link[stats]{median}}/\code{\link[stats]{mad}} in
  order to reduce computation time.
}
\examples{
## generate data
set.seed(1234)     # for reproducibility
x <- rnorm(10)     # standard normal
x[1] <- x[1] * 10  # introduce outlier

## standardize data
x
standardize(x)     # mean and sd
robStandardize(x)  # median and MAD
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{scale}}, \code{\link{sweep}}
}
\keyword{array}

