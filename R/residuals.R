# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Extract residuals from a sequence of regression models
#' 
#' Extract residuals from a sequence of regression models such as submodels 
#' along a robust or groupwise least angle regression sequence.
#' 
#' @method residuals seqModel
#' @aliases residuals.rlars residuals.grplars residuals.tslarsP
#' 
#' @param object  the model fit from which to extract residuals.
#' @param p  an integer giving the lag length for which to extract residuals 
#' (the default is to use the optimal lag length).
#' @param s  an integer vector giving the steps of the submodels for which to 
#' extract the residuals (the default is to use the optimal submodel).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one step.
#' @param \dots  for the \code{"tslars"} method, additional arguments to be 
#' passed down to the \code{"seqModel"} or \code{"optSeqModel"} method.  For 
#' the other methods, additional arguments are currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested residuals.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{residuals}}, \code{\link{rlars}}, 
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}}, 
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}}
#' 
#' @example inst/doc/examples/example-residuals.rlars.R
#' 
#' @keywords regression
#' 
#' @export

residuals.seqModel <- function(object, s = NA, drop = !is.null(s), ...) {
    getComponent(object, "residuals", s=s, drop=drop, ...)
}


#' @rdname residuals.seqModel
#' @method residuals optSeqModel
#' @export

residuals.optSeqModel <- function(object, ...) object$residuals


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
#' @seealso \code{\link[stats]{residuals}}, \code{\link{sparseLTS}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-residuals.sparseLTS.R
#' 
#' @keywords regression
#' 
#' @export

residuals.sparseLTS <- function(object, fit = c("reweighted", "raw", "both"), 
        standardized = FALSE, ...) {
    # initializations
    fit <- match.arg(fit)
    standardized <- isTRUE(standardized)
    # extract residuals and standardize if requested
    if(fit == "reweighted") {
        if(standardized) (object$residuals - object$center) / object$scale
        else object$residuals
    } else if(fit == "raw") {
        if(standardized) 
            (object$raw.residuals - object$raw.center) / object$raw.scale
        else object$raw.residuals
    } else {
        if(standardized) {
            reweighted <- (object$residuals - object$center) / object$scale
            raw <- (object$raw.residuals - object$raw.center) / object$raw.scale
            cbind(reweighted=reweighted, raw=raw)
        } else cbind(reweighted=object$residuals, raw=object$raw.residuals)
    }
}


#' @rdname residuals.sparseLTS
#' @method residuals sparseLTSGrid
#' @export

residuals.sparseLTSGrid <- function(object, s = NA, 
        fit = c("reweighted", "raw", "both"), 
        standardized = FALSE, drop = !is.null(s), ...) {
    ## extract residuals
    residuals <- getComponent(object, "residuals", s=s, fit=fit, drop=drop, ...)
    ## if requested, standardize residuals
    if(isTRUE(standardized)) {
        # extract center and scale estimates
        if(fit == "reweighted") {
            center <- object$center
            scale <- object$scale
        } else if(fit == "raw") {
            center <- object$raw.center
            scale <- object$raw.scale
        } else {
            center <- c(object$center, object$raw.center)
            scale <- c(object$scale, object$raw.scale)
        }
        # if requested, extract center and scale estimates for selected steps
        if(!is.null(s)) {
            center <- center[s]
            scale <- scale[s]
        }
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


#' @rdname residuals.sparseLTS
#' @method residuals optSparseLTSGrid
#' @export

residuals.optSparseLTSGrid <- residuals.sparseLTS
