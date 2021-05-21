# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @export
setupDiagnosticPlot <- function(object, ...) UseMethod("setupDiagnosticPlot")

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
    out <- list(data = data, facets = . ~ fit, qqLine = qql, q = q)
    class(out) <- "setupDiagnosticPlot"
  } else if (length(s) > 1) {
    ## extract the data from each requested step
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
    out <- list(data = data, facets = if(fit == "both") step ~ fit else ~step,
                qqLine = qql, q = q)
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
    out <- list(data = data, facets = . ~ fit, qqLine = qql, q = q)
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
  residuals <- residuals(object, s = s, fit = fit, standardized = TRUE)
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
  labs <- c("potential outlier", "regular observation")
  class <- ifelse(wt == 0, labs[1], labs[2])
  class <- factor(class, levels = labs)
  ## construct data frame containing main information
  data <- data.frame(index = seq_len(n), fitted = fitted, residual = residuals,
                     theoretical = theoretical, qqd = qqd, rd = rd, xyd = xyd,
                     weight = wt, leverage = leverage, diagnostics = class)
  ## construct return object
  out <- list(data = data, qqLine = as.data.frame(qql),
              q = data.frame(q = max(q, 2.5)))
  class(out) <- "setupDiagnosticPlot"
  out
}
