# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Extract coefficients from a sequence of regression models
#' 
#' Extract coefficients from a sequence of regression models, such as submodels 
#' along a robust least angle regression sequence, or sparse least trimmed 
#' squares regression models for a grid of values for the penalty parameter.
#' 
#' @method coef seqModel
#' @aliases coef.rlars
#' 
#' @param object  the model fit from which to extract coefficients.
#' @param s  for the \code{"seqModel"} method, an integer vector giving 
#' the steps of the submodels for which to extract coefficients (the default 
#' is to use the optimal submodel).  For the \code{"sparseLTS"} method, an 
#' integer vector giving the indices of the models for which to extract 
#' coefficients.  If \code{fit} is \code{"both"}, this can be a list with two 
#' components, with the first component giving the indices of the reweighted 
#' fits and the second the indices of the raw fits.  The default is to use the 
#' optimal model for each of the requested estimators.  Note that the optimal 
#' models may not correspond to the same value of the penalty parameter for the 
#' reweighted and the raw estimator.
#' @param fit  a character string specifying which coefficients to extract.  
#' Possible values are \code{"reweighted"} (the default) for the coefficients 
#' from the reweighted estimator, \code{"raw"} for the coefficients from the 
#' raw estimator, or \code{"both"} for the coefficients from both estimators.
#' @param zeros  a logical indicating whether to keep zero coefficients 
#' (\code{TRUE}, the default) or to omit them (\code{FALSE}).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one submodel.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested regression coefficients.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{coef}}, \code{\link{rlars}}, 
#' \code{\link{sparseLTS}}
#' 
#' @example inst/doc/examples/example-coef.R
#' 
#' @keywords regression
#' 
#' @export

coef.seqModel <- function(object, s = NA, zeros = TRUE, 
                          drop = !is.null(s), ...) {
  ## extract coefficients
  coef <- getComponent(object, "coefficients", s=s, drop=drop, ...)
  ## if requested, omit zero coefficients
  if(!isTRUE(zeros)) {
    if(is.null(dim(coef))) coef <- coef[coef != 0]
    else {
      keep <- apply(coef != 0, 1, any)
      coef <- coef[keep, , drop=FALSE]
    }
  }
  ## return coefficients
  coef
}


#' @rdname coef.seqModel
#' @method coef sparseLTS
#' @export

coef.sparseLTS <- function(object, s = NA, 
                           fit = c("reweighted", "raw", "both"), 
                           zeros = TRUE, drop = !is.null(s), ...) {
  ## extract coefficients
  coef <- getComponent(object, "coefficients", s=s, fit=fit, drop=drop, ...)
  ## if requested, omit zero coefficients
  if(!isTRUE(zeros)) {
    if(is.null(dim(coef))) coef <- coef[coef != 0]
    else {
      keep <- apply(coef != 0, 1, any)
      coef <- coef[keep, , drop=FALSE]
    }
  }
  ## return coefficients
  coef
}
