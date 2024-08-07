# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @noRd
weights.lmrob <- function(object, ...) object$rweights

#' @noRd
weights.lts <- function(object, ...) object$raw.weights

#' @noRd
weights.rlm <- function(object, ...) object$w


#' Extract outlier weights from sparse LTS regression models
#'
#' Extract binary weights that indicate outliers from sparse least trimmed
#' squares regression models.
#'
#' @param object  the model fit from which to extract outlier weights.
#' @param type  the type of weights to be returned.  Currently only robustness
#' weights are implemented (\code{"robustness"}).
#' @param s  an integer vector giving the indices of the models for which to
#' extract outlier weights.  If \code{fit} is \code{"both"}, this can be a list
#' with two components, with the first component giving the indices of the
#' reweighted fits and the second the indices of the raw fits.  The default is
#' to use the optimal model for each of the requested estimators.  Note that
#' the optimal models may not correspond to the same value of the penalty
#' parameter for the reweighted and the raw estimator.
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
#' @method weights sparseLTS
#' @export

weights.sparseLTS <- function(object, type = "robustness", s = NA,
                              fit = c("reweighted", "raw", "both"),
                              drop = !is.null(s), ...) {
  getComponent(object, "wt", s=s, fit=fit, drop=drop, ...)
}
