# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

weights.lmrob <- function(object, ...) object$rweights

weights.lts <- function(object, ...) object$raw.weights

weights.rlm <- function(object, ...) object$w

#' Extract robustness weights
#'
#' Extract robustness weights from robust penalized regression models.  For
#' sparse least trimmed squares regression models, those are binary weights
#' that indicate outliers.
#'
#' @param object  the model fit from which to extract robustness weights.
#' @param type  the type of weights to be returned.  Currently only robustness
#' weights are implemented (\code{"robustness"}).
#' @param s  an integer vector giving the indices of the models for which to
#' extract robustness weights.  For sparse least trimmed squares models, this
#' can be a list with two components (if \code{fit} is \code{"both"}), with the
#' first component giving the indices of the reweighted fits and the second the
#' indices of the raw fits.  The default is to use the optimal model for each
#' of the requested estimators.  Note that the optimal models may not
#' correspond to the same value of the penalty parameter for the reweighted and
#' the raw estimator.
#' @param fit  a character string specifying for which estimator to extract
#' outlier weights.  Possible values are \code{"reweighted"} (the default) for
#' weights indicating outliers from the reweighted fit, \code{"raw"} for
#' weights indicating outliers from the raw fit, or \code{"both"} for the
#' outlier weights from both estimators.
#' @param drop  a logical indicating whether to reduce the dimension to a
#' vector in case of only one model.
#' @param \dots  currently ignored.
#'
#' @return
#' A numeric vector or matrix containing the requested outlier weights.
#'
#' @note The weights are \eqn{1} for observations with reasonably small
#' residuals and \eqn{0} for observations with large residuals.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-weights.R
#'
#' @keywords regression
#'
#' @method weights penModel
#' @export

weights.penModel <- function(object, type = "robustness", s = NA,
                             drop = !is.null(s), ...) {
  getComponent(object, "rweights", s=s, drop=drop, ...)
}


#' @rdname weights.penModel
#' @method weights sparseLTS
#' @export

weights.sparseLTS <- function(object, type = "robustness", s = NA,
                              fit = c("reweighted", "raw", "both"),
                              drop = !is.null(s), ...) {
  getComponent(object, "wt", s=s, fit=fit, drop=drop, ...)
}
