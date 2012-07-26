\name{sparseLTSGrid}
\alias{print.sparseLTSGrid}
\alias{sparseLTSGrid}
\alias{sparseLTSGrid.default}
\alias{sparseLTSGrid.formula}
\title{Sparse least trimmed squares regression}
\usage{
  sparseLTSGrid(x, ...)

  \method{sparseLTSGrid}{formula} (formula, data, ...)

  \method{sparseLTSGrid}{default} (x, y, lambda,
    mode = c("lambda", "fraction"), ...,
    crit = c("BIC", "PE"), splits = foldControl(),
    cost = rtmspe, costArgs = list(),
    selectBest = c("hastie", "min"), seFactor = 1,
    seed = NULL, model = TRUE)
}
\arguments{
  \item{formula}{a formula describing the model.}

  \item{data}{an optional data frame, list or environment
  (or object coercible to a data frame by
  \code{\link{as.data.frame}}) containing the variables in
  the model.  If not found in data, the variables are taken
  from \code{environment(formula)}, typically the
  environment from which \code{sparseLTS} is called.}

  \item{x}{a numeric matrix containing the predictor
  variables.}

  \item{y}{a numeric vector containing the response
  variable.}

  \item{lambda}{a numeric vector of non-negative numeric
  values to be used as penalty parameter.}

  \item{mode}{a character string specifying the type of
  penalty parameter.  If \code{"lambda"}, \code{lambda}
  gives the grid of values for the penalty parameter
  directly.  If \code{"fraction"}, the smallest value of
  the penalty parameter that sets all coefficients to 0 is
  first estimated based on bivariate winsorization, then
  \code{lambda} gives the fractions of that estimate to be
  used (hence all values of \code{lambda} should be in the
  interval [0,1] in that case).}

  \item{crit}{a character string specifying the optimality
  criterion to be used for selecting the final model.
  Currently, only \code{"BIC"} for the Bayes information
  criterion is implemented.}

  \item{\dots}{additional arguments to be passed to
  \code{\link{sparseLTS}}.}

  \item{model}{a logical indicating whether the data
  \code{x} and \code{y} should be added to the return
  object.  If \code{intercept} is \code{TRUE}, a column of
  ones is added to \code{x} to account for the intercept.}
}
\value{
  An object of class \code{"sparseLTSGrid"} (inheriting
  from class \code{"seqModel"}) with the following
  components:

  \item{best}{an integer matrix in which each column
  contains the best subset of \eqn{h} observations found
  and used for computing the raw estimates with the
  corresponding penalty parameter.}

  \item{objective}{a numeric vector giving the values of
  the sparse LTS objective function, i.e., the
  \eqn{L_{1}}{L1} penalized sum of the \eqn{h} smallest
  squared residuals from the raw fits.}

  \item{coefficients}{a numeric matrix in which each column
  contains the coefficient estimates of the corresponding
  reweighted fit (including the intercept if
  \code{intercept} is \code{TRUE}).}

  \item{fitted.values}{a numeric matrix in which each
  column contains the fitted values of the response of the
  corresponding reweighted fit.}

  \item{residuals}{a numeric matrix in which each column
  contains the residuals of the response of the
  corresponding reweighted fit.}

  \item{center}{a numeric vector giving the robust center
  estimates of the residuals from the reweighted fits.}

  \item{scale}{a numeric vector giving the robust scale
  estimates of the residuals from the reweighted fits.}

  \item{lambda}{a numeric vector giving the values of the
  penalty parameter.}

  \item{intercept}{a logical indicating whether the model
  includes a constant term.}

  \item{alpha}{a numeric value giving the percentage of the
  residuals for which the \eqn{L_{1}}{L1} penalized sum of
  squares was minimized.}

  \item{quan}{the number \eqn{h} of observations used to
  compute the raw estimates.}

  \item{cnp2}{a numeric vector giving the consistency
  factors applied to the scale estimates of the residuals
  from the reweighted fits.}

  \item{wt}{an integer matrix in which each column contains
  binary weights that indicate outliers from the
  corresponding reweighted fit, i.e., the weights are
  \eqn{1} for observations with reasonably small reweighted
  residuals and \eqn{0} for observations with large
  reweighted residuals.}

  \item{df}{an integer vector giving the degrees of freedom
  of the obtained reweighted model fits, i.e., the number
  of nonzero coefficient estimates.}

  \item{raw.coefficients}{a numeric matrix in which each
  column contains the coefficient estimates of the
  corresponding raw fit (including the intercept if
  \code{intercept} is \code{TRUE}).}

  \item{raw.fitted.values}{a numeric matrix in which each
  column contains the fitted values of the response of the
  corresponding raw fit.}

  \item{raw.residuals}{a numeric matrix in which each
  column contains the residuals of the corresponding raw
  fit.}

  \item{raw.center}{a numeric vector giving the robust
  center estimates of the residuals from the raw fits.}

  \item{raw.scale}{a numeric vector giving the robust scale
  estimates of the residuals from the raw fits.}

  \item{raw.cnp2}{a numeric vector giving the consistency
  factors applied to the scale estimates of the residuals
  from the raw fits.}

  \item{raw.wt}{an integer matrix in which each column
  contains binary weights that indicate outliers of the
  corresponding raw fit, i.e., the weights used for the
  reweighted fits.}

  \item{crit}{a character string specifying the optimality
  criterion used for selecting the optimal model.}

  \item{critValues}{a numeric vector containing the values
  of the optimality criterion from the reweighted fits.}

  \item{sOpt}{an integer giving the optimal reweighted
  fit.}

  \item{raw.critValues}{a numeric vector containing the
  values of the optimality criterion from the raw fits.}

  \item{raw.sOpt}{an integer giving the optimal raw fit.}

  \item{x}{the predictor matrix (if \code{model} is
  \code{TRUE}).}

  \item{y}{the response variable (if \code{model} is
  \code{TRUE}).}

  \item{call}{the matched function call.}
}
\description{
  Compute least trimmed squares regression with an
  \eqn{L_{1}}{L1} penalty on the regression coefficients,
  which allows for sparse model estimates, over a grid of
  values for the penalty parameter.
}
\note{
  Package \pkg{robustHD} has a built-in back end for sparse
  least trimmed squares using the C++ library Armadillo.
  Another back end is available through package
  \pkg{sparseLTSEigen}, which uses the C++ library Eigen.
  The latter is faster, but not available on all platforms.
  For instance, \pkg{sparseLTSEigen} currently does not
  work on 32-bit \R for Windows.  In addition, there is
  currently no binary package for OS X available on CRAN
  due to problems with the PowerPC architecture.
  Nevertheless, OS X users with Intel machines can install
  \pkg{RcppEigen} and \pkg{sparseLTSEigen} from source if
  the standard \R developer tools are installed.
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

## fit sparse LTS models over a grid of values for lambda
frac <- seq(0.25, 0.05, by = -0.05)
sparseLTSGrid(x, y, lambda = frac, mode = "fraction")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{sparseLTS}},
  \code{\link[=coef.sparseLTS]{coef.sparseLTSGrid}},
  \code{\link[=fitted.sparseLTS]{fitted.sparseLTSGrid}},
  \code{\link[=plot.seqModel]{plot.sparseLTSGrid}},
  \code{\link{diagnosticPlot}},
  \code{\link[=predict.sparseLTS]{predict.sparseLTSGrid}},
  \code{\link[=residuals.sparseLTS]{residuals.sparseLTSGrid}},
  \code{\link[=wt.sparseLTS]{wt.sparseLTSGrid}},
}
\keyword{multivariate}
\keyword{robust}

