# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Extract residuals from sparse LTS regression models
#' 
#' Extract residuals from sparse least trimmed squares regression models.
#' 
#' @method residuals sparseLTS
#' 
#' @param object  the model fit from which to extract residuals.
#' @param s  an integer vector giving the indices of the models for which to 
#' extract residuals.  If \code{fit} is \code{"both"}, this can be a list 
#' with two components, with the first component giving the indices of the 
#' reweighted fits and the second the indices of the raw fits.  The default is 
#' to use the optimal model for each of the requested estimators.  Note that 
#' the optimal models may not correspond to the same value of the penalty 
#' parameter for the reweighted and the raw estimator.
#' @param fit  a character string specifying which residuals to extract.  
#' Possible values are \code{"reweighted"} (the default) for the residuals 
#' from the reweighted estimator, \code{"raw"} for the residuals from the raw 
#' estimator, or \code{"both"} for the residuals from both estimators.
#' @param standardized  a logical indicating whether the residuals should be 
#' standardized (the default is \code{FALSE}).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one model.
#' @param \dots  currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested (standardized) residuals.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{residuals}}, \code{\link{sparseLTS}}
#' 
#' @example inst/doc/examples/example-residuals.sparseLTS.R
#' 
#' @keywords regression
#' 
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
