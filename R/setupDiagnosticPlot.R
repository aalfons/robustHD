# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' Set up a diagnostic plot for a sequence of regression models
#'
#' Extract the fitted values and residuals of a sequence of regression models
#' (such as robust least angle regression models or sparse least trimmed
#' squares regression models) and other useful information for diagnostic
#' plots.
#'
#' @aliases setupDiagnosticPlot.rlars setupDiagnosticPlot.grplars
#' setupDiagnosticPlot.tslarsP
#'
#' @param object  the model fit from which to extract information.
#' @param s  for the \code{"seqModel"} method, an integer vector giving the
#' steps of the submodels from which to extract information (the default is to
#' use the optimal submodel).  For the \code{"sparseLTS"} method, an integer
#' vector giving the indices of the models from which to extract information
#' (the default is to use the optimal model for each of the requested fits).
#' @param p  an integer giving the lag length for which to extract information
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying from which fit to extract
#' information.  Possible values are \code{"reweighted"} (the default) to
#' convert the reweighted fit, \code{"raw"} to convert the raw fit, or
#' \code{"both"} to convert both fits.
#' @param covArgs  a list of arguments to be passed to
#' \code{\link[robustbase]{covMcd}} for computing robust Mahalanobis distances.
#' @param \dots  additional arguments to be passed to
#' \code{\link[robustbase]{covMcd}} can be specified directly instead of via
#' \code{covArgs}.
#'
#' @return  An object of class \code{"setupDiagnosticPlot"} with the following
#' components:
#' \describe{
#'   \item{\code{data}}{a data frame containing the columns listed below.
#'   \describe{
#'     \item{\code{step}}{the steps (for the \code{"seqModel"} method) or
#'     indices (for the \code{"sparseLTS"} method) of the models (only returned
#'     if more than one model is requested).}
#'     \item{\code{fit}}{the model fits (only returned if both the reweighted
#'     and raw fit are requested in the \code{"sparseLTS"} method).}
#'     \item{\code{index}}{the indices of the observations.}
#'     \item{\code{fitted}}{the fitted values.}
#'     \item{\code{residual}}{the standardized residuals.}
#'     \item{\code{theoretical}}{the corresponding theoretical quantiles from
#'     the standard normal distribution.}
#'     \item{\code{qqd}}{the absolute distances from a reference line through
#'     the first and third sample and theoretical quartiles.}
#'     \item{\code{rd}}{the robust Mahalanobis distances computed via the MCD
#'     (see \code{\link[robustbase]{covMcd}}).}
#'     \item{\code{xyd}}{the pairwise maxima of the absolute values of the
#'     standardized residuals and the robust Mahalanobis distances, divided by
#'     the respective other outlier detection cutoff point.}
#'     \item{\code{weight}}{the weights indicating regression outliers.}
#'     \item{\code{leverage}}{logicals indicating leverage points (i.e.,
#'     outliers in the predictor space).}
#'     \item{\code{Diagnostics}}{a factor with levels \code{"Potential outlier"}
#'     (potential regression outliers) and \code{"Regular observation"} (data
#'     points following the model).}
#'   }}
#'   \item{\code{qqLine}}{a data frame containing the intercepts and slopes of
#'   the respective reference lines to be displayed in residual Q-Q plots.}
#'   \item{\code{q}}{a data frame containing the quantiles of the Mahalanobis
#'   distribution used as cutoff points for detecting leverage points.}
#'   \item{\code{facets}}{default faceting formula for the diagnostic plots
#'   (only returned where applicable).}
#' }
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{diagnosticPlot}}, \code{\link{rlars}},
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
#' \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-setupDiagnosticPlot.R
#'
#' @keywords utilities
#'
#' @export

setupDiagnosticPlot <- function(object, ...) UseMethod("setupDiagnosticPlot")


#' @rdname setupDiagnosticPlot
#' @export
#'

setupDiagnosticPlot.seqModel <- function(object, s = NA, covArgs = list(...),
                                         ...) {
  ## initializations
  if (!object$robust) {
    stop("diagnostic plots not yet implemented for nonrobust methods")
  }
  # check the scale estimate
  scale <- getComponent(object, "scale", s = s)
  if (any(scale <= 0)) stop("residual scale equal to 0")
  # check if model data is available to compute robust MCD distances
  terms <- delete.response(object$terms)  # extract terms for model matrix
  x <- object$x
  if (is.null(x)) {
    x <- try(model.matrix(terms), silent = TRUE)
    if (inherits(x, "try-error")) {
      x <- NULL
      warning("model data not available")
    }
  }
  if (!is.null(x)) x <- removeIntercept(x)
  ## construct data frame with all information for plotting
  steps <- object$s
  if (length(steps) > 1) {
    ## check steps
    if (is.null(s)) s <- steps
    else if (isTRUE(is.na(s))) s <- getSOpt(object)  # defaults to optimal step
    else s <- checkSteps(s, sMin = steps[1], sMax = steps[length(steps)])
  } else s <- NA
  ## extract data for the requested steps
  if (length(s) > 1) {
    ## extract the information from each requested step
    out_list <- lapply(s, setupDiagnosticPlotSeqModelStep, object = object,
                       x = x, covArgs = covArgs)
    data_list <- lapply(out_list, "[[", "data")
    qql_list <- lapply(out_list, "[[", "qqLine")
    q_list <- lapply(out_list, "[[", "q")
    ## combine the information from the steps
    data <- cbind(step = rep.int(s, sapply(data_list, nrow)),
                  do.call(rbind, data_list))
    qql <- cbind(step = rep.int(s, sapply(qql_list, nrow)),
                 do.call(rbind, qql_list))
    q <- cbind(step = rep.int(s, sapply(q_list, nrow)),
               do.call(rbind, q_list))
    ## construct return object
    out <- list(data = data, qqLine = qql, q = q, facets = ~step)
    class(out) <- "setupDiagnosticPlot"
  } else {
    # extract the information from the selected step
    out <- setupDiagnosticPlotSeqModelStep(s, object = object, x = x,
                                           covArgs = covArgs)
  }
  ## return data frame
  out
}


## workhorse function for a single step from a sequence of regression models
setupDiagnosticPlotSeqModelStep <- function(s, object, x = NULL,
                                            covArgs = list()) {
  ## extract fitted values
  fitted <- fitted(object, s = s)
  ## extract standardized residuals
  residuals <- rstandard(object, s = s)
  n <- length(residuals)  # number of observations
  ## compute 0/1 outlier weights
  wt <- as.integer(abs(residuals) <= qnorm(0.9875))
  ## compute theoretical quantiles and distances from Q-Q reference line
  theoretical <- qqNorm(residuals)
  qql <- qqLine(residuals)  # Q-Q reference line
  qqd <- abs(residuals - qql$intercept - qql$slope * theoretical)
  ## compute MCD distances using selected variables
  # extract predictor matrix
  ok <- (is.na(s) || s > 0) && !is.null(x)
  if (ok) {
    # extract coefficients
    coefficients <- removeIntercept(coef(object, s = s))
    selected <- which(coefficients != 0)
    p <- length(selected)
    if (p == 0) {
      ok <- FALSE
      warning("all coefficients equal to 0")
    }
  }
  if (ok) {
    # compute distances
    rd <- try({
      xs <- x[, selected, drop = FALSE]
      callCovFun <- getCallFun(covArgs)
      mcd <- callCovFun(xs, fun = covMcd, args = covArgs)
      sqrt(mahalanobis(xs, mcd$center, mcd$cov))
    }, silent = TRUE)
    if(inherits(rd, "try-error")) {
      ok <- FALSE
      warning("robust distances cannot be computed")
    }
  }
  if (!ok) rd <- rep.int(NA, n)
  # take maximum of the distances in the x- and y-space, divided by the
  # respective other cutoff point
  q <- sqrt(qchisq(0.975, p))
  xyd <- pmax.int(abs(rd/2.5), abs(residuals/q))
  ## construct indicator variables for leverage points
  leverage <- rd > q
  ## classify data points
  labs <- c("Potential outlier", "Regular observation")
  class <- ifelse(wt == 0, labs[1], labs[2])
  class <- factor(class, levels = labs)
  ## construct data frame containing main information
  data <- data.frame(index = seq_len(n), fitted = fitted, residual = residuals,
                     theoretical = theoretical, qqd = qqd, rd = rd, xyd = xyd,
                     weight = wt, leverage = leverage, Diagnostics = class)
  ## construct return object
  out <- list(data = data, qqLine = as.data.frame(qql),
              q = data.frame(q = max(q, 2.5)))
  class(out) <- "setupDiagnosticPlot"
  out
}


#' @rdname setupDiagnosticPlot
#' @export

setupDiagnosticPlot.perrySeqModel <- function(object, ...) {
  # call method for component 'finalModel'
  setupDiagnosticPlot(object$finalModel, ...)
}


#' @rdname setupDiagnosticPlot
#' @export

setupDiagnosticPlot.tslars <- function(object, p, ...) {
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
  setupDiagnosticPlot(object$pFit[[p]], ...)
}


#' @rdname setupDiagnosticPlot
#' @export

setupDiagnosticPlot.sparseLTS <- function(object, s = NA,
                                          fit = c("reweighted", "raw", "both"),
                                          covArgs = list(...), ...) {
  ## initializations
  fit <- match.arg(fit)
  # check the scale estimate
  scale <- getComponent(object, "scale", s = s, fit = fit)
  if (any(scale <= 0)) stop("residual scale equal to 0")
  # check if model data is available to compute robust MCD distances
  terms <- delete.response(object$terms)  # extract terms for model matrix
  x <- object$x
  if (is.null(x)) {
    x <- try(model.matrix(terms), silent = TRUE)
    if(inherits(x, "try-error")) {
      x <- NULL
      warning("model data not available")
    }
  }
  if (!is.null(x) && object$intercept) x <- removeIntercept(x)
  ## construct data frame with all information for plotting
  lambda <- object$lambda
  bothOpt <- FALSE
  if (length(lambda) > 1) {
    ## check steps
    steps <- seq_along(lambda)
    if (is.null(s)) s <- steps
    else if (isTRUE(is.na(s))) {
      s <- getSOpt(object, fit = fit)  # defaults to optimal step
      if(fit == "both") {
        s <- unique(s)
        bothOpt <- length(s) == 2
      }
    } else s <- checkSteps(s, sMin = 1, sMax = length(steps))
  } else s <- NA
  ## extract data for the requested steps
  if (bothOpt) {
    # extract the data from the respecitve optimal lambda
    fits <- c("reweighted", "raw")
    ## recursive call for each fit
    reweighted <- setupDiagnosticPlotSparseLTSFit(object, s = s[1],
                                                  fit = "reweighted",
                                                  x = x)
    raw <- setupDiagnosticPlotSparseLTSFit(object, s = s[2], fit = "raw",
                                           x = x)
    ## combine data for Q-Q reference line
    qql <- data.frame(fit = factor(fits, levels = fits),
                      rbind(reweighted$qqLine, raw$qqLine),
                      row.names = NULL)
    ## combine data for cutoff chi-squared quantile
    q <- data.frame(fit = factor(fits, levels = fits),
                    rbind(reweighted$q, raw$q),
                    row.names = NULL)
    ## combine main information
    n <- c(nrow(reweighted$data), nrow(raw$data))
    data <- data.frame(fit = rep.int(factor(fits, levels = fits), n),
                       rbind(reweighted$data, raw$data),
                       row.names = NULL)
    ## construct return object
    out <- list(data = data, qqLine = qql, q = q, facets = . ~ fit)
    class(out) <- "setupDiagnosticPlot"
  } else if (length(s) > 1) {
    ## extract the information from each requested step
    out_list <- lapply(s, setupDiagnosticPlotSparseLTSStep, object = object,
                       fit = fit, x = x, covArgs = covArgs)
    data_list <- lapply(out_list, "[[", "data")
    qql_list <- lapply(out_list, "[[", "qqLine")
    q_list <- lapply(out_list, "[[", "q")
    ## combine the information from the steps
    data <- cbind(step = rep.int(s, sapply(data_list, nrow)),
                  do.call(rbind, data_list))
    qql <- cbind(step = rep.int(s, sapply(qql_list, nrow)),
                 do.call(rbind, qql_list))
    q <- cbind(step = rep.int(s, sapply(q_list, nrow)),
               do.call(rbind, q_list))
    ## construct return object
    facets <- if (fit == "both") step ~ fit else ~step
    out <- list(data = data, qqLine = qql, q = q, facets = facets)
    class(out) <- "setupDiagnosticPlot"
  } else {
    # extract the information from the selected step
    out <- setupDiagnosticPlotSparseLTSStep(s, object = object, fit = fit,
                                            x = x, covArgs = covArgs)
  }
  ## return object
  out
}


## workhorse functions

# setup diagnostic plots for a single sparse LTS step
setupDiagnosticPlotSparseLTSStep <- function(s, object, fit = "reweighted",
                                             x = NULL, covArgs = list()) {
  ## construct data frame with all information for plotting
  if (fit == "both") {
    fits <- c("reweighted", "raw")
    ## call workhorse function for each fit
    reweighted <- setupDiagnosticPlotSparseLTSFit(object, s = s,
                                                  fit = "reweighted",
                                                  x = x, covArgs = covArgs)
    raw <- setupDiagnosticPlotSparseLTSFit(object, s = s, fit = "raw",
                                           x = x, covArgs = covArgs)
    ## combine data for Q-Q reference line
    qql <- data.frame(fit = factor(fits, levels = fits),
                      rbind(reweighted$qqLine, raw$qqLine),
                      row.names = NULL)
    ## combine data for cutoff chi-squared quantile
    q <- data.frame(fit = factor(fits, levels = fits),
                    rbind(reweighted$q, raw$q),
                    row.names = NULL)
    ## combine main information
    n <- c(nrow(reweighted$data), nrow(raw$data))
    data <- data.frame(fit = rep.int(factor(fits, levels = fits), n),
                       rbind(reweighted$data, raw$data), row.names = NULL)
    ## construct return object
    out <- list(data = data, qqLine = qql, q = q, facets = . ~ fit)
    class(out) <- "setupDiagnosticPlot"
  } else {
    out <- setupDiagnosticPlotSparseLTSFit(object, s = s, fit = fit, x = x,
                                           covArgs = covArgs)
  }
  ## return object
  out
}

# setup diagnostic plots for a single sparse LTS fit
setupDiagnosticPlotSparseLTSFit <- function(object, s, fit = "reweighted",
                                            x = NULL, covArgs = list()) {
  ## extract fitted values
  fitted <- fitted(object, s = s, fit = fit)
  ## extract standardized residuals
  residuals <- rstandard(object, s = s, fit = fit)
  n <- length(residuals)  # number of observations
  ## extract outlier weights
  wt <- weights(object, type = "robustness", s = s, fit = fit)
  ## compute theoretical quantiles and distances from Q-Q reference line
  theoretical <- qqNorm(residuals)
  qql <- qqLine(residuals)  # Q-Q reference line
  qqd <- abs(residuals - qql$intercept - qql$slope * theoretical)
  ## compute MCD distances using selected variables
  # extract predictor matrix
  ok <- !is.null(x)
  if (ok) {
    # extract coefficients
    coefficients <- coef(object, s = s, fit = fit)
    if (object$intercept) coefficients <- removeIntercept(coefficients)
    selected <- which(coefficients != 0)
    p <- length(selected)
    if (p == 0) {
      ok <- FALSE
      warning("all coefficients equal to 0")
    }
  }
  if (ok) {
    # adjust alpha since MCD computes subset size depending on n and p
    h <- object$quan
    n2 <- (n + p + 1) %/% 2
    alpha <- pmin((h - 2*n2 + n) / (2 * (n - n2)), 1)
    # check fraction for subset size
    if (alpha < 0.5) {
      alpha <- 0.5
      warning(sprintf("cannot compute MCD with h = %d; using h = %d",
                      object$quan, h.alpha.n(alpha, n, p)))
    }
    # compute distances
    rd <- try({
      xs <- x[, selected, drop = FALSE]
      covArgs$alpha <- alpha
      callCovFun <- getCallFun(covArgs)
      mcd <- callCovFun(xs, fun = covMcd, args = covArgs)
      if (fit == "reweighted") {
        center <- mcd$center
        cov <- mcd$cov
      } else {
        center <- mcd$raw.center
        cov <- mcd$raw.cov
      }
      sqrt(mahalanobis(xs, center, cov))
    }, silent = TRUE)
    if(inherits(rd, "try-error")) {
      ok <- FALSE
      warning("robust distances cannot be computed")
    }
  }
  if (!ok) rd <- rep.int(NA_real_, n)
  # take maximum of the distances in the x- and y-space, divided by the
  # respective other cutoff point
  q <- sqrt(qchisq(0.975, p))
  xyd <- pmax.int(abs(rd/2.5), abs(residuals/q))
  ## construct indicator variables for leverage points
  leverage <- rd > q
  ## classify data points
  labs <- c("Potential outlier", "Regular observation")
  class <- ifelse(wt == 0, labs[1], labs[2])
  class <- factor(class, levels = labs)
  ## construct data frame containing main information
  data <- data.frame(index = seq_len(n), fitted = fitted, residual = residuals,
                     theoretical = theoretical, qqd = qqd, rd = rd, xyd = xyd,
                     weight = wt, leverage = leverage, Diagnostics = class)
  ## construct return object
  out <- list(data = data, qqLine = as.data.frame(qql),
              q = data.frame(q = max(q, 2.5)))
  class(out) <- "setupDiagnosticPlot"
  out
}


#' @rdname setupDiagnosticPlot
#' @export

setupDiagnosticPlot.perrySparseLTS <- function(object, ...) {
  # call method for component 'finalModel'
  setupDiagnosticPlot(object$finalModel, ...)
}
