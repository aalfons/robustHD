# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @export
setupCoefPlot <- function(object, ...) UseMethod("setupCoefPlot")

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

#' @export
setupCoefPlot.tslars <- function(object, p, ...) {
  ## check lag length
  if (missing(p) || !is.numeric(p) || length(p) == 0) p <- object$pOpt
  if (length(p) > 1) {
    warning("multiple lag lengths not yet supported")
    p <- p[1]
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

#' @export
setupCoefPlot.sparseLTS <- function(object,
                                    fit = c("reweighted", "raw", "both"),
                                    zeros = FALSE, labels = NULL, ...) {
  # initializations
  fit <- match.arg(fit)
  zeros <- isTRUE(zeros)
  coef <- removeIntercept(t(coef(object, s = NULL, fit = fit)))
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
      out$facets = . ~ fit
    } else {
      out$labels <- data.frame(lambda = rep.int(unname(lambda[sMax]), m),
                               step = rep.int(steps[sMax], m),
                               df = rep.int(unname(dfSMax), m),
                               coefficient = unname(coef[sMax, ]),
                               label = labels)
    }
  }
  # add class and return object
  class(out) <- "setupCoefPlot"
  out
}
