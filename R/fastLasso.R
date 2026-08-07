# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Lasso regression on a subset
#'
#' Wrapper function for a barebones implementation in \proglang{C++} to
#' compute the lasso on a subset of the data for a fixed value of the penalty
#' parameter via the LARS algorithm. It is mainly used internally by
#' \code{\link{sparseLTS}()} to compute the reweighted fit, but also for the
#' raw fit in the rare case of no trimming.
#'
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a non-negative numeric value to be used as penalty parameter.
#' @param subset  an integer vector giving the subset on which to compute the
#' lasso.  If \code{NULL}, the lasso is computed on the full data set.
#' @param normalize  a logical indicating whether the predictor variables
#' should be normalized to have unit \eqn{L_{2}}{L2} norm (the default is
#' \code{TRUE}).  Note that normalization is performed on the subsample
#' rather than the full data set.
#' @param intercept  a logical indicating whether a constant term should be
#' included in the model (the default is \code{TRUE}).
#' @param eps  a small positive numeric value used to determine whether the
#' variability within a variable is too small (an effective zero).
#' @param use.Gram  a logical indicating whether the Gram matrix of the
#' explanatory variables should be precomputed on the subsample (the default is
#' \code{TRUE}).  If the number of variables is large, computation may be
#' faster when this is set to \code{FALSE}.
#' @param drop a logical indicating whether relevant components of the output
#' should be returned as vectors (\code{TRUE}, the default) or matrices
#' (\code{FALSE}).
#' @param raw  argument for internal use that is not expected to be used by
#' others. (If \code{TRUE}, somewhat different output is returned that is
#' relevant for the raw fit without trimming in \code{\link{sparseLTS}()}.)
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{coefficients}}{a numeric vector or matrix containing the
#'   coefficient estimates (as determined by argument \code{drop}).}
#'   \item{\code{fitted.values}}{a numeric vector or matrix containing the
#'   fitted values (as determined by argument \code{drop}).}
#'   \item{\code{residuals}}{a numeric vector or matrix containing the
#'   residuals (as determined by argument \code{drop}).}
#' }
#'
#' @note
#' As of version 0.8.5, this function is exported for use by other developers
#' due to popular demand. Use this function only if you know what you are
#' doing! There is no error handling, and the underlying \proglang{C++} code
#' may crash the \proglang{R} session if improper arguments are supplied.
#'
#' @seealso \code{\link{sparseLTS}}
#'
#' @keywords regression
#'
#' @author Andreas Alfons
#'
#' @references
#' Efron, B., Hastie, T., Johnstone, I. and Tibshirani, R. (2003) Least angle
#' regression. \emph{The Annals of Statistics}, \bold{32}(2), 407--499.
#' \doi{10.1214/009053604000000067}
#'
#' @export

fastLasso <- function(x, y, lambda, subset = NULL, normalize = TRUE,
                      intercept = TRUE, eps = .Machine$double.eps,
                      use.Gram = TRUE, drop = TRUE, raw = FALSE) {
  # initializations
  intercept <- isTRUE(intercept)
  use.Gram <- isTRUE(use.Gram)
  drop <- isTRUE(drop)
  raw <- isTRUE(raw)
  # compute lasso
  if(raw) {
    # call C++ function
    fit <- .Call("R_testLasso", R_x=x, R_y=y, R_lambda=lambda,
                 R_initial=seq_along(y), R_normalize=normalize,
                 R_intercept=intercept, R_eps=eps, R_useGram=use.Gram,
                 PACKAGE = "robustHD")

    # prepare object for raw lasso fit
    coef <- fit$coefficients
    res <- fit$residuals
    if(drop) {
      # drop the dimension of the components
      coef <- drop(coef)
      res <- drop(res)
    }
    center <- mean(res)
    scale <- sqrt(mean((res-center)^2))
    fit <- list(best=fit$indices, coefficients=coef, residuals=res,
                objective=fit$crit, center=center, scale=scale)
  } else {
    # check subset
    if(is.null(subset)) {
      useSubset <- FALSE
      subset <- integer()
    } else useSubset <- TRUE
    # call C++ function
    fit <- .Call("R_fastLasso", R_x=x, R_y=y, R_lambda=lambda,
                 R_useSubset=useSubset, R_subset=subset,
                 R_normalize=normalize, R_intercept=intercept,
                 R_eps=eps, R_useGram=use.Gram,
                 PACKAGE = "robustHD")
    if(drop) fit <- lapply(fit, drop)  # drop the dimension of the components
  }
  # return lasso fit
  fit
}
