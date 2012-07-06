# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## class 'rlm'

weights.rlm <- function(object, ...) object$w


#' Extract outlier weights from sparse LTS regression models
#' 
#' Extract binary weights that indicate outliers from sparse least trimmed 
#' squares regression models.
#' 
#' @param object  the model fit from which to extract outlier weights.
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
#' @seealso \code{\link{sparseLTS}}, \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-wt.sparseLTS.R
#' 
#' @keywords regression
#' 
#' @export

wt <- function(object, ...) UseMethod("wt")


#' @rdname wt
#' @method wt sparseLTS
#' @export

wt.sparseLTS <- function(object, fit = c("reweighted", "raw", "both"), ...) {
    fit <- match.arg(fit)
    switch(fit, reweighted=object$wt, raw=object$raw.wt,
        both=cbind(reweighted=object$wt, raw=object$raw.wt))
}


#' @rdname wt
#' @method wt sparseLTSGrid
#' @export

wt.sparseLTSGrid <- function(object, s, fit = c("reweighted", "raw", "both"), 
        drop = !is.null(s), ...) {
    ## initializations
    fit <- match.arg(fit)
    ## extract weights
    if(fit == "reweighted") wt <- object$wt
    else if(fit == "raw") wt <- object$raw.wt
    else {
        wt <- list(reweighted=object$wt, raw=object$raw.wt)
        wt <- mapply(function(x, n) {
                colnames(x) <- paste(n, colnames(x), sep=".")
                x
            }, wt, names(wt), SIMPLIFY=FALSE)
        wt <- do.call(cbind, wt)
    }
    ## check selected steps and extract corresponding weights
    sMax <- length(object$lambda)
    if(missing(s)) 
        s <- switch(fit, reweighted=object$sOpt, raw=object$raw.sOpt, 
            both=c(reweighted=object$sOpt, raw=sMax+object$raw.sOpt))
    else if(!is.null(s)) {
        if(fit == "both" && is.list(s)) {
            s <- rep(s, length.out=2)
            s <- lapply(s, checkSteps, sMin=1, sMax=sMax)
            s <- c(s[[1]], sMax+s[[2]])
        } else {
            s <- checkSteps(s, sMin=1, sMax=sMax)
            if(fit == "both") s <- c(s, sMax+s)
        }
    }
    if(!is.null(s)) wt <- wt[, s, drop=FALSE]  # selected steps
    ## return weights
    if(isTRUE(drop)) drop(wt) else wt
}
