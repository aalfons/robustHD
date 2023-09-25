# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Extract residuals from a sequence of regression models
#'
#' Extract residuals from a sequence of regression models, such as submodels
#' along a robust or groupwise least angle regression sequence, or sparse least
#' trimmed squares regression models for a grid of values for the penalty
#' parameter.
#'
#' @method residuals seqModel
#' @aliases residuals.rlars residuals.grplars residuals.tslarsP
#'
#' @param object  the model fit from which to extract residuals.
#' @param p  an integer giving the lag length for which to extract residuals
#' (the default is to use the optimal lag length).
#' @param s  for the \code{"seqModel"} method, an integer vector giving the
#' steps of the submodels for which to extract the residuals (the default is to
#' use the optimal submodel).  For the \code{"sparseLTS"} method, an integer
#' vector giving the indices of the models for which to extract residuals.  If
#' \code{fit} is \code{"both"}, this can be a list with two components, with
#' the first component giving the indices of the reweighted fits and the second
#' the indices of the raw fits.  The default is to use the optimal model for
#' each of the requested estimators.  Note that the optimal models may not
#' correspond to the same value of the penalty parameter for the reweighted and
#' the raw estimator.
#' @param fit  a character string specifying which residuals to extract.
#' Possible values are \code{"reweighted"} (the default) for the residuals
#' from the reweighted estimator, \code{"raw"} for the residuals from the raw
#' estimator, or \code{"both"} for the residuals from both estimators.
#' @param standardized  a logical indicating whether the residuals should be
#' standardized (the default is \code{FALSE}). Note that this argument is
#' \bold{deprecated} and may be removed as soon as the next version. Use
#' \code{\link[=rstandard.seqModel]{rstandard}} instead to extract standardized
#' residuals.
#' @param drop  a logical indicating whether to reduce the dimension to a
#' vector in case of only one step.
#' @param \dots  for the \code{"tslars"} method, additional arguments to be
#' passed down to the \code{"seqModel"} method.  For the other methods,
#' additional arguments are currently ignored.
#'
#' @return
#' A numeric vector or matrix containing the requested residuals.
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link[stats]{residuals}}, \code{\link[=rstandard.seqModel]{rstandard}}
#'
#' \code{\link{rlars}}, \code{\link{grplars}}, \code{\link{rgrplars}},
#' \code{\link{tslarsP}}, \code{\link{rtslarsP}}, \code{\link{tslars}},
#' \code{\link{rtslars}}, \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-residuals.R
#'
#' @keywords regression
#'
#' @import stats
#' @export

residuals.seqModel <- function(object, s = NA, standardized = FALSE,
                               drop = !is.null(s), ...) {
  # check for deprecated argument 'standardized'
  if(isTRUE(standardized)) {
    warning("Argument 'standardized' is deprecated.\n",
            "Use the 'rstandard' method instead to extract standardized ",
            "residuals.")
    rstandard(object, s = s, drop = drop, ...)
  } else {
    # extract residuals
    residuals <- getComponent(object, "residuals", s = s, drop = FALSE, ...)
    # drop dimension if requested and return residuals
    if(isTRUE(drop)) dropCol(residuals) else residuals
  }
}


#' @rdname residuals.seqModel
#' @method residuals tslars
#' @export

residuals.tslars <- function(object, p, ...) {
  ## initializations
  # check lag length
  if(missing(p) || !is.numeric(p) || length(p) == 0) {
    p <- object$pOpt
  } else p <- p[1]
  pMax <- object$pMax
  if(p < 1) {
    p <- 1
    warning("lag length too small, using lag length 1")
  } else if(p > pMax) {
    p <- pMax
    warning(sprintf("lag length too large, using maximum lag length %d", p))
  }
  ## extract residuals for specified lag length
  residuals(object$pFit[[p]], ...)
}


#' @rdname residuals.seqModel
#' @method residuals perrySeqModel
#' @export
residuals.perrySeqModel <- function(object, ...) {
  finalModel <- object$finalModel
  if(is.null(finalModel)) stop("final model not available")
  residuals(finalModel, ...)
}


#' @rdname residuals.seqModel
#' @method residuals sparseLTS
#' @export

residuals.sparseLTS <- function(object, s = NA,
                                fit = c("reweighted", "raw", "both"),
                                standardized = FALSE, drop = !is.null(s),
                                ...) {
  # check for deprecated argument 'standardized'
  if(isTRUE(standardized)) {
    warning("Argument 'standardized' is deprecated.\n",
            "Use the 'rstandard' method instead to extract standardized ",
            "residuals.")
    rstandard(object, s = s, fit = fit, drop = drop, ...)
  } else {
    # extract residuals
    getComponent(object, "residuals", s = s, fit = fit, drop = drop, ...)
  }
}
