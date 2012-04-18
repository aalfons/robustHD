# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Extract coefficients from a sequence of regression models
#' 
#' Extract coefficients from a sequence of regression models such as submodels 
#' along a robust or groupwise least angle regression sequence.
#' 
#' @method coef seqModel
#' @aliases coef.rlars coef.grplars coef.tslarsP
#' 
#' @param object  the model fit from which to extract coefficients.
#' @param p  an integer giving the lag length for which to extract coefficients 
#' (the default is to use the optimal lag length).
#' @param s  an integer vector giving the steps of the submodels for which to 
#' extract coefficients (the default is to use the optimal submodel).
#' @param zeros  a logical indicating whether to keep zero coefficients 
#' (\code{TRUE}, the default) or to drop them (\code{FALSE}).
#' @param \dots  for the \code{"tslars"} method, additional arguments to be 
#' passed down to the \code{"seqModel"} method.  For the \code{"seqModel"} 
#' method, additional arguments are currently ignored.
#' 
#' @return  
#' If only one submodel is requested, a numeric vector containing the 
#' corresponding regression coefficients.
#' 
#' If multiple submodels are requested, a numeric matrix in which each column 
#' contains the regression coefficients of the corresponding submodel.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{coef}}, \code{\link{rlars}}, 
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}}, 
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}}
#' 
#' @example inst/doc/examples/example-coef.rlars.R
#' 
#' @keywords regression
#' 
#' @export

coef.seqModel <- function(object, s, zeros = TRUE, ...) {
    ## extract coefficients
    coef <- getComponent(object, "coefficients", s=s, ...)
    ## if requested, drop zero coefficients in case of a single step
    if(!isTRUE(zeros)) {
        if(is.null(dim(coef))) {
            coef <- coef[coef != 0]
        } else {
            keep <- apply(coef != 0, 1, any)
            coef <- coef[keep, , drop=FALSE]
        }
    }
    ## return coefficients
    coef
}


#' @rdname coef.seqModel
#' @method coef tslars
#' @export

coef.tslars <- function(object, p, ...) {
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
    ## extract coefficients for specified lag length
    coef(object$pFit[[p]], ...)
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
#' (\code{TRUE}, the default) or to drop them (\code{FALSE}).
#' @param \dots  currently ignored.
#' 
#' @return  
#' If coefficients for only one model are requested, they are returned in the 
#' form of a numeric vector.
#' 
#' Otherwise a numeric matrix is returned in which each column contains the 
#' coefficients of the corresponding model.
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
    fit <- match.arg(fit)
    coef <- switch(fit,
        reweighted=object$coefficients,
        raw=object$raw.coefficients,
        both=cbind(reweighted=object$coefficients, raw=object$raw.coefficients))
    if(!isTRUE(zeros)) {
        # drop zero coefficients
        if(is.null(dim(coef))) {
            coef <- coef[coef != 0]
        } else {
            keep <- apply(coef != 0, 1, any)
            coef <- coef[keep, , drop=FALSE]
        }
    }
    coef
}


#' @rdname coef.sparseLTS
#' @method coef sparseLTSGrid
#' @export

coef.sparseLTSGrid <- function(object, s, fit = c("reweighted", "raw", "both"), 
        zeros = TRUE, ...) {
    ## initializations
    fit <- match.arg(fit)
    ## extract coefficients
    if(fit == "reweighted") {
        coef <- object$coefficients
    } else if(fit == "raw") {
        coef <- object$raw.coefficients
    } else {
        coef <- list(reweighted=object$coefficients, raw=object$raw.coefficients)
        coef <- mapply(function(x, n) {
                colnames(x) <- paste(n, colnames(x), sep=".")
                x
            }, coef, names(coef), SIMPLIFY=FALSE)
        coef <- do.call(cbind, coef)
    }
    ## check selected steps and extract corresponding coefficients
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
    if(!is.null(s)) coef <- coef[, s]  # coefficients for selected steps
    ## if requested, drop zero coefficients in case of a single step
    if(!isTRUE(zeros)) {
        if(is.null(dim(coef))) {
            coef <- coef[coef != 0]
        } else {
            keep <- apply(coef != 0, 1, any)
            coef <- coef[keep, , drop=FALSE]
        }
    }
    ## return coefficients
    coef
}
