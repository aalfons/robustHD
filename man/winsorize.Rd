% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/winsorize.R
\name{winsorize}
\alias{winsorize}
\alias{winsorize.default}
\alias{winsorize.matrix}
\alias{winsorize.data.frame}
\title{Data cleaning by winsorization}
\usage{
winsorize(x, ...)

\method{winsorize}{default}(
  x,
  standardized = FALSE,
  centerFun = median,
  scaleFun = mad,
  const = 2,
  return = c("data", "weights"),
  ...
)

\method{winsorize}{matrix}(
  x,
  standardized = FALSE,
  centerFun = median,
  scaleFun = mad,
  const = 2,
  prob = 0.95,
  tol = .Machine$double.eps^0.5,
  return = c("data", "weights"),
  ...
)

\method{winsorize}{data.frame}(x, ...)
}
\arguments{
\item{x}{a numeric vector, matrix or data frame to be cleaned.}

\item{\dots}{for the generic function, additional arguments to be passed
down to methods.  For the \code{"data.frame"} method, additional arguments
to be passed down to the \code{"matrix"} method.  For the other methods,
additional arguments to be passed down to
\code{\link[=standardize]{robStandardize}}.}

\item{standardized}{a logical indicating whether the data are already
robustly standardized.}

\item{centerFun}{a function to compute a robust estimate for the center to
be used for robust standardization (defaults to
\code{\link[stats]{median}}).  Ignored if \code{standardized} is \code{TRUE}.}

\item{scaleFun}{a function to compute a robust estimate for the scale to
be used for robust standardization (defaults to \code{\link[stats]{mad}}).
Ignored if \code{standardized} is \code{TRUE}.}

\item{const}{numeric; tuning constant to be used in univariate
winsorization (defaults to 2).}

\item{return}{character string; if \code{standardized} is \code{TRUE},
this specifies the type of return value.  Possible values are \code{"data"}
for returning the cleaned data, or \code{"weights"} for returning data
cleaning weights.}

\item{prob}{numeric; probability for the quantile of the
\eqn{\chi^{2}}{chi-squared} distribution to be used in multivariate
winsorization (defaults to 0.95).}

\item{tol}{a small positive numeric value used to determine singularity
issues in the computation of correlation estimates based on bivariate
winsorization (see \code{\link{corHuber}}).}
}
\value{
If \code{standardize} is \code{TRUE} and \code{return} is \code{"weights"},
a set of data cleaning weights.  Multiplying each observation of the
standardized data by the corresponding weight yields the cleaned
standardized data.

Otherwise an object of the same type as the original data \code{x}
containing the cleaned data is returned.
}
\description{
Clean data by means of winsorization, i.e., by shrinking outlying
observations to the border of the main part of the data.
}
\details{
The borders of the main part of the data are defined on the scale of the
robustly standardized data.  In the univariate case, the borders are given
by \eqn{+/-}\code{const}, thus a symmetric distribution is assumed.  In the
multivariate case, a normal distribution is assumed and the data are
shrunken towards the boundary of a tolerance ellipse with coverage
probability \code{prob}.  The boundary of this ellipse is thereby given by
all points that have a squared Mahalanobis distance equal to the quantile of
the \eqn{\chi^{2}}{chi-squared} distribution given by \code{prob}.
}
\note{
Data cleaning weights are only meaningful for standardized data.  In
the general case, the data need to be standardized first, then the data
cleaning weights can be computed and applied to the standardized data, after
which the cleaned standardized data need to be backtransformed to the
original scale.
}
\examples{
## generate data
set.seed(1234)     # for reproducibility
x <- rnorm(10)     # standard normal
x[1] <- x[1] * 10  # introduce outlier

## winsorize data
x
winsorize(x)

}
\references{
Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust linear model
selection based on least angle regression. \emph{Journal of the American
Statistical Association}, \bold{102}(480), 1289--1299.
\doi{10.1198/016214507000000950}
}
\seealso{
\code{\link{corHuber}}
}
\author{
Andreas Alfons, based on code by Jafar A. Khan, Stefan Van Aelst and
Ruben H. Zamar
}
\keyword{robust}
