\name{getScale}
\alias{getScale}
\alias{getScale.seqModel}
\alias{getScale.sparseLTS}
\title{Extract the residual scale of a robust regression model}
\usage{
  getScale(x, ...)

  \method{getScale}{seqModel} (x, s = NA, ...)

  \method{getScale}{sparseLTS} (x, s = NA,
    fit = c("reweighted", "raw", "both"), ...)
}
\arguments{
  \item{x}{the model fit from which to extract the robust
  residual scale estimate.}

  \item{s}{for the \code{"seqModel"} method, an integer
  vector giving the steps of the submodels for which to
  extract the robust residual scale estimate (the default
  is to use the optimal submodel).  For the
  \code{"sparseLTS"} method, an integer vector giving the
  indices of the models from which to extract the robust
  residual scale estimate.  If \code{fit} is \code{"both"},
  this can be a list with two components, with the first
  component giving the indices of the reweighted fits and
  the second the indices of the raw fits.  The default is
  to use the optimal model for each of the requested
  estimators.  Note that the optimal models may not
  correspond to the same value of the penalty parameter for
  the reweighted and the raw estimator.}

  \item{fit}{a character string specifying from which fit
  to extract the robust residual scale estimate.  Possible
  values are \code{"reweighted"} (the default) for the
  residual scale of the reweighted fit, \code{"raw"} for
  the residual scale of the raw fit, or \code{"both"} for
  the residual scale of both fits.}

  \item{\dots}{additional arguments to be passed down to
  methods.}
}
\value{
  A numeric vector or matrix giving the robust residual
  scale estimates for the requested model fits.
}
\description{
  Extract the robust scale estimate of the residuals from a
  robust regression model.
}
\details{
  Methods are implemented for models of class
  \code{"lmrob"} (see \code{\link[robustbase]{lmrob}}),
  \code{"lts"} (see \code{\link[robustbase]{ltsReg}}),
  \code{"rlm"} (see \code{\link[MASS]{rlm}}),
  \code{"seqModel"} (see \code{\link{rlars}}) and
  \code{"sparseLTS"} (see \code{\link{sparseLTS}}).  The
  default method computes the MAD of the residuals.
}
\examples{
data("coleman")
fit <- lmrob(Y ~ ., data=coleman)
getScale(fit)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[=AIC.seqModel]{AIC}},
  \code{\link[robustbase]{lmrob}},
  \code{\link[robustbase]{ltsReg}},
  \code{\link[MASS]{rlm}}, \code{\link{rlars}},
  \code{\link{sparseLTS}}
}
\keyword{regression}

