# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Extract fitted values from a sequence of regression models
#' 
#' Extract fitted values from a sequence of regression models such as submodels 
#' along a robust or groupwise least angle regression sequence.
#' 
#' @method fitted seqModel
#' @aliases fitted.rlars fitted.grplars fitted.tslarsP
#' 
#' @param object  the model fit from which to extract fitted values.
#' @param p  an integer giving the lag length for which to extract fitted 
#' values (the default is to use the optimal lag length).
#' @param s  an integer vector giving the steps of the submodels for which to 
#' extract the fitted values (the default is to use the optimal submodel).
#' @param drop  a logical indicating whether to reduce the dimension to a 
#' vector in case of only one step.
#' @param \dots  for the \code{"tslars"} method, additional arguments to be 
#' passed down to the \code{"seqModel"} method.  For the \code{"seqModel"} 
#' method, additional arguments are currently ignored.
#' 
#' @return  
#' A numeric vector or matrix containing the requested fitted values.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{fitted}}, \code{\link{rlars}}, 
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}}, 
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}}
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


#' @rdname fitted.seqModel
#' @method fitted tslars
#' @export

fitted.tslars <- function(object, p, ...) {
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
    ## extract fitted values for specified lag length
    fitted(object$pFit[[p]], ...)
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
    switch(fit,
        reweighted=object$fitted.values,
        raw=object$raw.fitted.values,
        both=cbind(reweighted=object$fitted.values, 
            raw=object$raw.fitted.values))
}


#' @rdname fitted.sparseLTS
#' @method fitted sparseLTSGrid
#' @export

fitted.sparseLTSGrid <- function(object, s, 
        fit = c("reweighted", "raw", "both"), 
        drop = !is.null(s), ...) {
    ## initializations
    fit <- match.arg(fit)
    ## extract fitted values
    if(fit == "reweighted") {
        fitted <- object$fitted.values
    } else if(fit == "raw") {
        fitted <- object$raw.fitted.values
    } else {
        fitted <- list(reweighted=object$fitted.values, 
            raw=object$raw.fitted.values)
        fitted <- mapply(function(x, n) {
                colnames(x) <- paste(n, colnames(x), sep=".")
                x
            }, fitted, names(fitted), SIMPLIFY=FALSE)
        fitted <- do.call(cbind, fitted)
    }
    ## check selected steps and extract corresponding fitted values
    sMax <- length(object$lambda)
    if(missing(s)) {
        s <- switch(fit, reweighted=object$sOpt, raw=object$raw.sOpt, 
            both=c(reweighted=object$sOpt, raw=sMax+object$raw.sOpt))
    } else if(!is.null(s)) {
        if(fit == "both" && is.list(s)) {
            s <- rep(s, length.out=2)
            s <- lapply(s, checkSteps, sMin=1, sMax=sMax)
            s <- c(s[[1]], sMax+s[[2]])
        } else {
            s <- checkSteps(s, sMin=1, sMax=sMax)
            if(fit == "both") s <- c(s, sMax+s)
        }
    }
    if(!is.null(s)) fitted <- fitted[, s, drop=FALSE]  # selected steps
    ## return fitted values
    if(isTRUE(drop)) drop(fitted) else fitted
}
