# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Extract residuals from a sequence of regression models
#' 
#' Extract residuals from a sequence of regression models, such as submodels 
#' along a robust least angle regression sequence, or sparse least trimmed 
#' squares regression models for a grid of values for the penalty parameter.
#' 
#' @method residuals seqModel
#' @aliases residuals.rlars
#' 
#' @param object  the model fit from which to extract residuals.
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
#' standardized (the default is \code{FALSE}).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one step.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested residuals.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{residuals}}, \code{\link{rlars}}, 
#' \code{\link{sparseLTS}}
#' 
#' @example inst/doc/examples/example-residuals.R
#' 
#' @keywords regression
#' 
#' @export

residuals.seqModel <- function(object, s = NA, standardized = FALSE, 
                               drop = !is.null(s), ...) {
  ## extract residuals
  residuals <- getComponent(object, "residuals", s=s, drop=drop, ...)
  ## if requested, standardize residuals
  if(isTRUE(standardized)) {
    if(object$robust) {
      # extract scale estimates
      scale <- getComponent(object, "scale", s=s, ...)
      # standardize selected residuals
      if(is.null(dim(residuals))) residuals <- residuals / scale
      else residuals <- sweep(residuals, 2, scale, "/", check.margin=FALSE)
    } else stop("not implemented yet")
  }
  ## return residuals
  residuals
}


#' @rdname residuals.seqModel
#' @method residuals sparseLTS
#' @export

residuals.sparseLTS <- function(object, s = NA, 
                                fit = c("reweighted", "raw", "both"), 
                                standardized = FALSE, drop = !is.null(s), 
                                ...) {
  ## extract residuals
  residuals <- getComponent(object, "residuals", s=s, fit=fit, drop=drop, ...)
  ## if requested, standardize residuals
  if(isTRUE(standardized)) {
    # extract center and scale estimates
    center <- getComponent(object, "center", s=s, fit=fit, ...)
    scale <- getComponent(object, "scale", s=s, fit=fit, ...)
    # standardize selected residuals
    if(is.null(dim(residuals))) residuals <- (residuals - center) / scale
    else {
      residuals <- x <- sweep(residuals, 2, center, check.margin=FALSE)
      residuals <- sweep(residuals, 2, scale, "/", check.margin=FALSE)
    }
  }
  ## return residuals
  residuals
}
