# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## get best subset
getBest <- function(x, ...) UseMethod("getBest")

getBest.sparseLTSGrid <- function(x, s = NULL, drop = !is.null(s), ...) {
    best <- x$best
    if(!is.null(s)) best <- best[, s, drop=FALSE]
    if(drop) dropCol(best) else best
}


## get residual scale
getScale <- function(x, ...) UseMethod("getScale")

getScale.sparseLTSGrid <- function(x, s = NA, fit = "reweighted", ...) {
    getVector(x, "scale", s=s, fit=fit, ...)
}

getScale.optSparseLTSGrid <- function(x, fit = "reweighted", ...) {
    switch(fit, reweighted=x$scale, raw=x$raw.scale, 
        both=c(reweighted=x$scale, raw=x$raw.scale))
}


## get optimal step
getSOpt <- function(x, ...) UseMethod("getSOpt")

getSOpt.sparseLTSGrid <- function(x, fit = "reweighted", ...) {
    switch(fit, reweighted=x$sOpt, raw=x$raw.sOpt, 
        both=c(reweighted=x$sOpt, raw=x$raw.sOpt))
}

getSOpt.optSparseLTSGrid <- function(x, fit = "reweighted", ...) {
    sOpt <- x$critValues$best
    if(fit != "both") sOpt <- sOpt[fit]
    sOpt
}


## get a vector component for certain steps of a model sequence

# generic function
getVector <- function(x, component, ...) UseMethod("getVector")

# method for class "seqModel"
getVector.seqModel <- function(x, component, s = NA, ...) {
    comp <- x[[component]]  # extract component
    if(!is.null(s)) {
        if(isTRUE(is.na(s))) s <- x$sOpt  # optimal step size as default
        comp <- comp[s+1]                 # extract selected steps
    }
    comp  # return component
}

# method for class "sparseLTSGrid"
getVector.sparseLTSGrid <- function(x, component, s = NA, 
        fit = "reweighted", ...) {
    # initializations
    if(fit != "reweighted") raw.component <- paste("raw", component, sep=".")
    # extract component
    comp <- switch(fit, reweighted=x[[component]], raw=x[[raw.component]], 
        both=unlist(list(reweighted=x[[component]], raw=x[[raw.component]])))
    # extract selected parts of the component if requested
    if(!is.null(s)) {
        sMax <- length(x$lambda)
        if(isTRUE(is.na(s))) {
            s <- getSOpt(x, fit=fit)  # optimal step size as default
            if(fit == "both") s[2] <- sMax + s[2]
        } else if(fit == "both") s <- c(s, sMax+s)
        comp <- comp[s]  # extract selected steps
    }
    # return component
    comp
}


## get a matrix component for certain steps of a model sequence

# generic function
getMatrix <- function(x, component, ...) UseMethod("getMatrix")

# method for class "seqModel"
getMatrix.seqModel <- function(x, component, s = NA, drop = !is.null(s), ...) {
    comp <- x[[component]]  # extract component
    if(!is.null(s)) {
        if(isTRUE(is.na(s))) s <- x$sOpt  # optimal step size as default
        comp <- comp[, s+1, drop=FALSE]   # extract selected steps
    }
    if(drop) dropCol(comp) else comp  # drop dimension if requested
}

# method for class "sparseLTSGrid"
getMatrix.sparseLTSGrid <- function(x, component, s = NA, 
        fit = "reweighted", drop = !is.null(s), ...) {
    # initializations
    if(fit != "reweighted") raw.component <- paste("raw", component, sep=".")
    # extract component
    comp <- switch(fit, reweighted=x[[component]], raw=x[[raw.component]], 
        both={
            reweighted <- x[[component]]
            colnames(reweighted) <- paste("reweighted", 
                colnames(reweighted), sep=".")
            raw <- x[[raw.component]]
            colnames(raw) <- paste("raw", colnames(raw), sep=".")
            cbind(reweighted, raw)
        })
    # extract selected parts of the component if requested
    if(!is.null(s)) {
        sMax <- length(x$lambda)
        if(isTRUE(is.na(s))) {
            s <- getSOpt(x, fit=fit)  # optimal step size as default
            if(fit == "both") s[2] <- sMax + s[2]
        } else if(fit == "both") s <- c(s, sMax+s)
        comp <- comp[, s, drop=FALSE]  # extract selected steps
    }
    # drop dimension if requested
    if(drop) dropCol(comp) else comp
}


## get a matrix component for certain steps of the model sequence
## this is used for accessors that are exported to the namespace 
## (coefficients, fitted values, residuals, ...), so checks for the 
## arguments are necessary

# generic function
getComponent <- function(x, component, ...) UseMethod("getComponent")

# method for class "seqModel"
getComponent.seqModel <- function(x, component, 
        s = NA, drop = !is.null(s), ...) {
    # extract component
    comp <- x[[component]]
    # check selected steps and extract corresponding parts of the component
    if(!is.null(s)) {
        if(isTRUE(is.na(s))) s <- x$sOpt  # defaults to optimal step size
        else s <- checkSteps(s, sMin=0, sMax=ncol(comp)-1)  # check steps
    }
    # extract selected steps
    # the extra check for NULL is necessary since component 'sOpt' is not 
    # available in internal object used in prediction error estimation
    if(!is.null(s)) comp <- comp[, s+1, drop=FALSE]
    # drop dimension if requested and return component
    if(isTRUE(drop)) dropCol(comp) else comp
}

# method for class "sparseLTSGrid"
getComponent.sparseLTSGrid <- function(x, component, s = NA, haveRaw = TRUE, 
        fit = c("reweighted", "raw", "both"), drop = !is.null(s), ...) {
    # initializations
    fit <- match.arg(fit)
    if(fit != "reweighted") raw.component <- paste("raw", component, sep=".")
    # extract component
    comp <- switch(fit, reweighted=x[[component]], raw=x[[raw.component]], 
        both={
            reweighted <- x[[component]]
            colnames(reweighted) <- paste("reweighted", 
                colnames(reweighted), sep=".")
            raw <- x[[raw.component]]
            colnames(raw) <- paste("raw", colnames(raw), sep=".")
            cbind(reweighted, raw)
        })
    # check selected steps and extract corresponding parts of the component
    if(!is.null(s)) {
        sMax <- length(x$lambda)
        if(isTRUE(is.na(s))) {
            s <- getSOpt(x, fit=fit)  # optimal step size as default
            if(fit == "both") s[2] <- sMax + s[2]
        } else if(fit == "both" && is.list(s)) {
            s <- rep(s, length.out=2)
            s <- lapply(s, checkSteps, sMin=1, sMax=sMax)
            s <- c(s[[1]], sMax+s[[2]])
        } else {
            s <- checkSteps(s, sMin=1, sMax=sMax)
            if(fit == "both") s <- c(s, sMax+s)
        }
        comp <- comp[, s, drop=FALSE]  # extract selected steps
    }
    # drop dimension if requested and return component
    if(isTRUE(drop)) dropCol(comp) else comp
}
