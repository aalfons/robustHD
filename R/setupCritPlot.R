# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @export
setupCritPlot <- function(object, ...) UseMethod("setupCritPlot")


#' @export
setupCritPlot.seqModel <- function(object, which = c("line", "dot"), ...) {
  ## extract optimality criterion
  crit <- object$crit
  if (is.null(crit)) stop("optimality criterion data not available")
  which <- match.arg(which)
  ## construct object containing relevant information for plotting
  setupCritPlot(crit, which = which, s = object$s, ...)
}


#' @export
setupCritPlot.tslars <- function(object, p, ...) {
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
  ## call plot function for specified lag length
  setupCritPlot(object$pFit[[p]], ...)
}


#' @export
setupCritPlot.sparseLTS <- function(object, which = c("line", "dot"),
                                    fit = c("reweighted", "raw", "both"),
                                    ...) {
  ## initializations
  crit <- object$crit
  if (is.null(crit)) stop("optimality criterion data not available")
  which <- match.arg(which)
  fit <- match.arg(fit)
  select <- if (fit == "both") NULL else fit
  ## construct object containing relevant information for plotting
  setupCritPlot(crit, which = which,
                tuning = data.frame(lambda = object$lambda),
                select = select)
}


## workhorse function to prepare BIC values for plotting
# object ... object of class "bicSelect"
# s ........ integer vector of step numbers for the sequence of models
# tuning ... data frame containing additional information such as values of a
#            tuning parameter
# select ... indicates columns of a BIC matrix to keep
setupCritPlot.bicSelect <- function(object, which = "line", s = NULL,
                                    tuning = NULL, select = NULL, ...) {
  # extract BIC values and make sure they are in matrix form
  bic <- object$values
  if (!is.null(dim(bic)) && !is.null(select)) bic <- bic[, select, drop = FALSE]
  bic <- as.data.frame(bic)
  # define index of the submodels
  d <- dim(bic)
  # check data frames for BIC and additional information
  if (any(d == 0)) stop("BIC data has no rows or columns")
  if (is.null(s)) s <- seq_len(d[1])
  else if (!is.numeric(s) || length(s) != d[1]) {
    stop(sprintf("'s' must be a vector of length %d", d[1]))
  }
  if (!is.null(tuning) && (!is.data.frame(tuning) || nrow(tuning) != d[1])) {
    stop(sprintf("'tuning' must be a data frame with %d rows", d[1]))
  }
  # construct data frame with information on optimality criterion
  if (!is.null(tuning) && ncol(tuning) == 1) fits <- tuning[, , drop = TRUE]
  else fits <- s
  if (d[2] == 1) {
    # combine information on the sequence with BIC values
    bic <- cbind(fits, bic)
    names(bic) <- c("Fit", "BIC")
    facets <- NULL
  } else {
    # reshape BIC matrix and add information on the sequence
    bic <- mapply(function(val, nam) {
      data.frame(Fit = fits, Name = rep.int(nam, d[1]), BIC = val)
    }, val = bic, nam = names(bic), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    bic <- do.call(rbind, bic)
    facets <- . ~ Name
  }
  # construct return object
  out <- list(data = bic, which = which, includeSE = FALSE,
              mapping = aes_string(x = "Fit", y = "BIC",
                                   ymin = "BIC", ymax = "BIC"))
  if (!is.null(facets)) out$facets <- facets
  if (!is.null(tuning)) out$tuning <- tuning
  class(out) <- c("setupBICPlot", "setupCritPlot")
  out
}


#' @export
setupCritPlot.perrySeqModel <- function(object,
                                        which = c("line", "dot",
                                                  "box", "density"),
                                        ...) {
  # initializations
  if (object$splits$R > 1) which <- match.arg(which)
  else {
    choices <- eval(formals()[["which"]])
    if (identical(which, choices)) which <- "line"
    else which <- match.arg(which, c("line", "dot"))
  }
  ## call setupPerryPlot() and add superclass
  out <- setupPerryPlot(object, which = which)
  class(out) <- c(class(out), "setupCritPlot")
  out
}


#' @export
setupCritPlot.perrySparseLTS <- function(object,
                                         which = c("line", "dot",
                                                   "box", "density"),
                                         fit = c("reweighted", "raw", "both"),
                                         ...) {
  ## initializations
  if (object$splits$R > 1) which <- match.arg(which)
  else {
    choices <- eval(formals()[["which"]])
    if (identical(which, choices)) which <- "line"
    else which <- match.arg(which, c("line", "dot"))
  }
  fit <- match.arg(fit)
  if (fit == "both") fit <- NULL
  ## call setupPerryPlot() and add superclass
  out <- setupPerryPlot(object, which = which, select = fit)
  ## if only one fit is to be plotted, the plot is facetted by default
  if (!is.null(fit)) out$facets <- NULL
  class(out) <- c(class(out), "setupCritPlot")
  out
}
