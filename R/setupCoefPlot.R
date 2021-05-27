# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up a coefficient plot of a sequence of regression models
#'
#' Extract the relevent information for a plot of the coefficients for a
#' sequence of regression models, such as submodels along a robust or groupwise
#' least angle regression sequence, or sparse least trimmed squares regression
#' models for a grid of values for the penalty parameter.
#'
#' @aliases setupCoefPlot.rlars setupCoefPlot.grplars
#' setupCoefPlot.tslarsP
#'
#' @param object  the model fit from which to extract information.
#' @param zeros  a logical indicating whether predictors that never enter the
#' model and thus have zero coefficients should be included in the plot
#' (\code{TRUE}) or omitted (\code{FALSE}, the default).  This is useful if the
#' number of predictors is much larger than the number of observations, in
#' which case many coefficients are never nonzero.
#' @param labels  an optional character vector containing labels for the
#' predictors.  Information on labels can be suppressed by setting this to
#' \code{NA}.
#' @param p  an integer giving the lag length for which to extract information
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying for which estimator to extract
#' information.  Possible values are \code{"reweighted"} (the default) for
#' the reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for
#' both estimators.
#' @param \dots  additional arguments to be passed down.
#'
#' @return  An object inheriting from class \code{"setupCoefPlot"} with the
#' following components:
#' \describe{
#'   \item{\code{coefficients}}{a data frame containing the following columns:
#'   \describe{
#'     \item{\code{fit}}{the model fit for which the coefficient is computed
#'     (only returned if both the reweighted and raw fit are requested in the
#'     \code{"sparseLTS"} method).}
#'     \item{\code{lambda}}{the value of the penalty parameter for which the
#'     coefficient is computed (only returned for the \code{"sparseLTS"}
#'     method).}
#'     \item{\code{step}}{the step along the sequence for which the coefficient
#'     is computed.}
#'     \item{\code{df}}{the degrees of freedom of the submodel along the
#'     sequence for which the coefficient is computed.}
#'     \item{\code{coefficient}}{the value of the coefficient.}
#'     \item{\code{variable}}{a character string specifying to which variable
#'     the coefficient belongs.}
#'   }
#'   }
#'   \item{\code{abscissa}}{a character string specifying available options for
#'   what to plot on the \eqn{x}-axis}
#'   \item{\code{lambda}}{a numeric vector giving the values of the penalty
#'   parameter. (only returned for the \code{"sparseLTS"} method).}
#'   \item{\code{step}}{an integer vector containing the steps for which
#'   submodels along the sequence have been computed.}
#'   \item{\code{df}}{an integer vector containing the degrees of freedom of
#'   the submodels along the sequence (i.e., the number of estimated
#'   coefficients; only returned for the \code{"seqModel"} method).}
#'   \item{\code{includeLabels}}{a logical indicating whether information on
#'   labels for the variables should be included in the plot.}
#'   \item{\code{labels}}{a data frame containing the following columns (not
#'   returned if information on labels is suppressed):
#'   \describe{
#'     \item{\code{fit}}{the model fit for which the coefficient is computed
#'     (only returned if both the reweighted and raw fit are requested in the
#'     \code{"sparseLTS"} method).}
#'     \item{\code{lambda}}{the smallest value of the penalty parameter
#'     (only returned for the \code{"sparseLTS"}
#'     method).}
#'     \item{\code{step}}{the last step along the sequence.}
#'     \item{\code{df}}{the degrees of freedom of the last submodel along the
#'     sequence.}
#'     \item{\code{coefficient}}{the value of the coefficient.}
#'     \item{\code{label}}{the label of the corresponding variable to be
#'     displayed in the plot.}
#'   }
#'   }
#'   \item{\code{facets}}{default faceting formula for the plots (only
#'   returned if both estimators are requested in the \code{"sparseLTS"}
#'   method).}
#' }
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{coefPlot}}, \code{\link{rlars}},
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
#' \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-setupCoefPlot.R
#'
#' @keywords utilities
#'
#' @export

setupCoefPlot <- function(object, ...) UseMethod("setupCoefPlot")


#' @rdname setupCoefPlot
#' @export

setupCoefPlot.seqModel <- function(object, zeros = FALSE, labels = NULL, ...) {
  # initializations
  zeros <- isTRUE(zeros)
  coef <- removeIntercept(t(coef(object, s = NULL)))
  m <- ncol(coef)  # number of variables
  # check labels
  includeLabels <- !isTRUE(is.na(labels))
  if (includeLabels) {
    if (is.null(labels)) labels <- defaultLabels(object)
    else {
      ok <- (is.numeric(labels) || is.character(labels)) && length(labels) == m
      if (!ok) {
        stop(sprintf("'labels' should be a character vector of length %d", m))
      }
    }
  }
  # prepare coefficients and labels
  sigmaX <- object$sigmaX
  if (!zeros) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop = FALSE]
    sigmaX <- sigmaX[keep]
    if (includeLabels) labels <- labels[keep]
  }
  # standardize coefficients
  coef <- sweep(coef, 2, sigmaX, "/", check.margin = FALSE)
  # prepare other information
  m <- ncol(coef)          # number of variables
  vn <- colnames(coef)     # variable names
  steps <- object$s        # step numbers
  nsteps <- length(steps)  # number of steps
  df <- object$df          # degrees of freedom
  # construct data frame for coefficients
  coefData <- data.frame(step = rep.int(steps, m),
                         df = rep.int(df, m), coefficient = as.numeric(coef),
                         variable = factor(rep(vn, each = nsteps), levels = vn))
  # construct return object
  out <- list(coefficients = coefData, abscissa = c("step", "df"),
              steps = steps, df = df, includeLabels = includeLabels)
  # add data frame for labels if requested
  if(includeLabels) {
    out$labels <- data.frame(step = rep.int(steps[nsteps], m),
                             df = rep.int(df[nsteps], m),
                             coefficient = coef[nsteps, ],
                             label = labels)
  }
  # add class and return object
  class(out) <- "setupCoefPlot"
  out
}


#' @rdname setupCoefPlot
#' @export

setupCoefPlot.tslars <- function(object, p, ...) {
  ## check lag length
  if (missing(p) || !is.numeric(p) || length(p) == 0) p <- object$pOpt
  if (length(p) > 1) {
    p <- p[1]
    warning(sprintf("multiple lag lengths not yet supported; using p = %d", p))
  }
  pMax <- object$pMax
  if (p < 1) {
    p <- 1
    warning("lag length too small, using lag length 1")
  } else if (p > pMax) {
    p <- pMax
    warning(sprintf("lag length too large, using maximum lag length %d", p))
  }
  ## call method for specified lag length
  setupCoefPlot(object$pFit[[p]], ...)
}


#' @rdname setupCoefPlot
#' @export

setupCoefPlot.sparseLTS <- function(object,
                                    fit = c("reweighted", "raw", "both"),
                                    zeros = FALSE, labels = NULL, ...) {
  # initializations
  fit <- match.arg(fit)
  zeros <- isTRUE(zeros)
  coef <- removeIntercept(t(coef(object, s = NULL, fit = fit)))
  m <- ncol(coef)  # number of variables
  df <- getComponent(object, "df", s = NULL, fit = fit)
  # check labels
  includeLabels <- !isTRUE(is.na(labels))
  if (includeLabels) {
    if (is.null(labels)) labels <- defaultLabels(object)
    else {
      ok <- (is.numeric(labels) || is.character(labels)) && length(labels) == m
      if (!ok) {
        stop(sprintf("'labels' should be a character vector of length %d", m))
      }
    }
  }
  # prepare coefficients and labels
  if (!zeros) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop = FALSE]
    if (includeLabels) labels <- labels[keep]
  }
  # check if predictor data is available to compute scale estimates
  x <- object$x
  if (is.null(x)) {
    x <- try(model.matrix(object$terms), silent = TRUE)
    if (inherits(x, "try-error"))
      stop("scale estimates of predictor variables not available")
  }
  x <- removeIntercept(x)
  if (!zeros) x <- x[, keep, drop = FALSE]
  #   obtain scale estimates for predictors
  n <- nrow(x)  # number of observations
  sigmaX <- apply(x, 2, function(x) {
    # standardize data
    xs <- robStandardize(x, fallback = TRUE)
    # detect good data points
    ok <- which(abs(xs) < qnorm(0.9875))
    nOk <- length(ok)
    # compute consistency factor
    if (nOk < n) {
      qn <- qnorm((nOk+n) / (2*n))  # quantile for consistency factor
      cdelta <- 1 / sqrt(1-(2*n)/(nOk/qn)*dnorm(qn))
    } else cdelta <- 1  # consistency factor not necessary
    # compute standard deviation of good data points and multiply with
    # consistency factor
    sd(x[ok]) * cdelta
  })
  # standardize coeffients
  coef <- sweep(coef, 2, sigmaX, "/", check.margin = FALSE)
  # prepare other information
  m <- ncol(coef)             # number of variables
  vn <- colnames(coef)        # variable names
  lambda <- object$lambda      # tuning parameters
  steps <- seq_along(lambda)  # step numbers
  sMax <- length(steps)       # number of steps
  # construct data frame for coefficients
  if(fit == "both") {
    fits <- c("reweighted", "raw")
    coefData <- data.frame(
      fit = rep.int(factor(rep(fits, each = sMax), levels = fits), m),
      lambda = rep.int(lambda, 2*m), step = rep.int(steps, 2*m),
      df = rep.int(df, m), coefficient = as.numeric(coef),
      variable = factor(rep(vn, each = 2*sMax), levels = vn)
    )
  } else {
    coefData <- data.frame(
      lambda = rep.int(lambda, m), step = rep.int(steps, m),
      df = rep.int(df, m), coefficient = as.numeric(coef),
      variable = factor(rep(vn, each = sMax), levels = vn)
    )
  }
  # construct return object
  out <- list(coefficients = coefData, abscissa = c("lambda", "step"),
              lambda = lambda, steps = steps, includeLabels = includeLabels)
  # add data frame for labels if requested
  if(includeLabels) {
    dfSMax <- getComponent(object, "df", s = sMax, fit = fit)
    if(fit == "both") {
      out$labels <- data.frame(fit = rep.int(factor(fits, levels = fits), m),
                               lambda = rep.int(unname(lambda[sMax]), 2*m),
                               step = rep.int(steps[sMax], 2*m),
                               df = rep.int(unname(dfSMax), m),
                               coefficient = unname(coef[sMax, ]),
                               label = rep(labels, each = 2))
    } else {
      out$labels <- data.frame(lambda = rep.int(unname(lambda[sMax]), m),
                               step = rep.int(steps[sMax], m),
                               df = rep.int(unname(dfSMax), m),
                               coefficient = unname(coef[sMax, ]),
                               label = labels)
    }
  }
  if (fit == "both")  out$facets = . ~ fit
  # add class and return object
  class(out) <- "setupCoefPlot"
  out
}
