# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Extract coefficients from a sequence of regression models
#' 
#' Extract coefficients from a sequence of regression models such as submodels 
#' along a robust least angle regression sequence.
#' 
#' @method coef seqModel
#' @aliases coef.rlars
#' 
#' @param object  the model fit from which to extract coefficients.
#' @param s  an integer vector giving the steps of the submodels for which to 
#' extract coefficients (the default is to use the optimal submodel).
#' @param zeros  a logical indicating whether to keep zero coefficients 
#' (\code{TRUE}, the default) or to omit them (\code{FALSE}).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one step.
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested regression coefficients.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{coef}}, \code{\link{rlars}}
#' 
#' @example inst/doc/examples/example-coef.rlars.R
#' 
#' @keywords regression
#' 
#' @export

coef.seqModel <- function(object, s, zeros = TRUE, drop = !is.null(s), ...) {
    ## extract coefficients
    if(missing(s) && missing(drop)) drop <- TRUE
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


#' Extract coefficients from sparse LTS regression models
#' 
#' Extract coefficients from sparse least trimmed squares regression models.  
#' 
#' @method coef sparseLTS
#' 
#' @param object  the model fit from which to extract coefficients.
#' @param s  an integer vector giving the indices of the models for which to 
#' extract coefficients.  If \code{fit} is \code{"both"}, this can be a list 
#' with two components, with the first component giving the indices of the 
#' reweighted fits and the second the indices of the raw fits.  The default is 
#' to use the optimal model for each of the requested estimators.  Note that 
#' the optimal models may not correspond to the same value of the penalty 
#' parameter for the reweighted and the raw estimator.
#' @param fit  a character string specifying which coefficients to extract.  
#' Possible values are \code{"reweighted"} (the default) for the coefficients 
#' from the reweighted estimator, \code{"raw"} for the coefficients from the 
#' raw estimator, or \code{"both"} for the coefficients from both estimators.
#' @param zeros  a logical indicating whether to keep zero coefficients 
#' (\code{TRUE}, the default) or to omit them (\code{FALSE}).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one model.
#' @param \dots  currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested regression coefficients.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{coef}}, \code{\link{sparseLTS}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-coef.sparseLTS.R
#' 
#' @keywords regression
#' 
#' @export

coef.sparseLTS <- function(object, fit = c("reweighted", "raw", "both"), 
        zeros = TRUE, ...) {
    ## extract coefficients
    fit <- match.arg(fit)
    coef <- switch(fit,
        reweighted=object$coefficients, raw=object$raw.coefficients,
        both=cbind(reweighted=object$coefficients, raw=object$raw.coefficients))
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


#' @rdname coef.sparseLTS
#' @method coef sparseLTSGrid
#' @export

coef.sparseLTSGrid <- function(object, s, 
        fit = c("reweighted", "raw", "both"), 
        zeros = TRUE, drop = !is.null(s), ...) {
    ## extract coefficients
    if(missing(s) && missing(drop)) drop <- TRUE
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
