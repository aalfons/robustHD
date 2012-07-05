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
#' @param \dots  additional arguments are currently ignored.
#' 
#' @return  
#' If only one submodel is requested, a numeric vector containing the 
#' corresponding fitted values.
#' 
#' If multiple submodels are requested, a numeric matrix in which each column 
#' contains the fitted values of the corresponding submodel.
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

fitted.seqModel <- function(object, s, drop = TRUE, ...) {
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
#' @param \dots  currently ignored.
#' 
#' @return  
#' If fitted values for only one model are requested, they are returned in the 
#' form of a numeric vector.
#' 
#' Otherwise a numeric matrix is returned in which each column contains the 
#' fitted values of the corresponding model.
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
        drop = TRUE, ...) {
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
    if(!is.null(s)) fitted <- fitted[, s, drop=drop]  # selected steps
    ## return fitted values
    fitted
}
