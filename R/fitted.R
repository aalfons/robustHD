# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Extract fitted values from a sequence of regression models
#' 
#' Extract fitted values from a sequence of regression models such as submodels 
#' along a robust least angle regression sequence.
#' 
#' @method fitted seqModel
#' @aliases fitted.rlars
#' 
#' @param object  the model fit from which to extract fitted values.
#' @param s  an integer vector giving the steps of the submodels for which to 
#' extract the fitted values (the default is to use the optimal submodel).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one step.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested fitted values.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{fitted}}, \code{\link{rlars}}
#' 
#' @example inst/doc/examples/example-fitted.rlars.R
#' 
#' @keywords regression
#' 
#' @export

fitted.seqModel <- function(object, s, drop = !is.null(s), ...) {
    if(missing(s) && missing(drop)) drop <- TRUE
    getComponent(object, "fitted.values", s=s, drop=drop, ...)
}


#' Extract fitted values from sparse LTS regression models
#' 
#' Extract fitted values from sparse least trimmed squares regression models.
#' 
#' @method fitted sparseLTS
#' 
#' @param object  the model fit from which to extract fitted values.
#' @param s  an integer vector giving the indices of the models for which to 
#' extract fitted values.  If \code{fit} is \code{"both"}, this can be a list 
#' with two components, with the first component giving the indices of the 
#' reweighted fits and the second the indices of the raw fits.  The default is 
#' to use the optimal model for each of the requested estimators.  Note that 
#' the optimal models may not correspond to the same value of the penalty 
#' parameter for the reweighted and the raw estimator.
#' @param fit  a character string specifying which fitted values to extract.  
#' Possible values are \code{"reweighted"} (the default) for the fitted values 
#' from the reweighted estimator, \code{"raw"} for the fitted values from the 
#' raw estimator, or \code{"both"} for the fitted values from both estimators.
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one model.
#' @param \dots  currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested fitted values.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats:fitted.values]{fitted}}, \code{\link{sparseLTS}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-fitted.sparseLTS.R
#' 
#' @keywords regression
#' 
#' @export

fitted.sparseLTS <- function(object, fit = c("reweighted", "raw", "both"), 
        ...) {
    fit <- match.arg(fit)
    switch(fit, reweighted=object$fitted.values, raw=object$raw.fitted.values,
        both=cbind(reweighted=object$fitted.values, 
            raw=object$raw.fitted.values))
}


#' @rdname fitted.sparseLTS
#' @method fitted sparseLTSGrid
#' @export

fitted.sparseLTSGrid <- function(object, s, 
        fit = c("reweighted", "raw", "both"), 
        drop = !is.null(s), ...) {
    if(missing(s) && missing(drop)) drop <- TRUE
    getComponent(object, "fitted.values", s=s, fit=fit, drop=drop, ...)
}
