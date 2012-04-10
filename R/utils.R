# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## add default column names to matrix
addColnames <- function(x) {
    # 'x' needs to be a matrix
    if(is.null(colnames(x))) colnames(x) <- paste("x", seq_len(ncol(x)), sep="")
    x
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
    if(!check || is.na(match("(Intercept)", colnames(x)))) {
        cbind("(Intercept)"=rep.int(1, nrow(x)), x)
    } else x
}

## remove intercept column from design matrix
#removeIntercept <- function(x, pos) {
#    if(missing(pos)) {
#        pos <- match("(Intercept)", colnames(x), nomatch = 0)
#        if(pos > 0) x[, -pos, drop=FALSE] else x
#    } else x[, -pos, drop=FALSE]
#}
removeIntercept <- function(x, pos) {
    haveVector <- is.null(dim(x))
    if(missing(pos)) {
        names <- if(haveVector) names(x) else colnames(x)
        pos <- match("(Intercept)", names, nomatch = 0)
    }
    if(pos > 0) {
        if(haveVector) x[-pos] else x[, -pos, drop=FALSE]
    } else x
}

### backtransform regression coefficients to original scale (including intercept)
#backtransform <- function(beta, muY, sigmaY, mu, sigma) {
#    apply(beta, 2, 
#        function(b) {
#            b <- b * sigmaY / sigma
#            a <- muY - sum(b * mu)  # intercept
#            c("(Intercept)"=a, b)
#        })
#}

## check maximum step for robust and groupwise LARS
# sequence variables as long as there are more observations than predictors
# fit models as long as there are twice as many observations than predictors
checkSMax <- function(sMax, n, p) {
    sMax <- rep(sMax, length.out=2)
    if(!isTRUE(is.finite(sMax[1])) || !isTRUE(sMax[1] <= min(p, n-1))) {
        sMax[1] <- min(p, n-1)
    }
#    if(!isTRUE(is.finite(sMax[2])) || !isTRUE(sMax[2] <= min(p, n-1))) {
#        sMax[2] <- min(p, n-1)
#    }
    if(!isTRUE(is.finite(sMax[2]))) {
        sMax[2] <- NA
    } else {
        if(!isTRUE(sMax[2] <= min(p, n-1))) sMax[2] <- min(p, n-1)
        if(sMax[2] > sMax[1]) sMax[2] <- sMax[1]
    }
    sMax
}

## check steps for coef(), fitted(), residuals(), predict(), ... methods
checkSteps <- function(s, sMin, sMax) {
    if(!is.numeric(s) || length(s) == 0 || any(s < sMin) || any(s > sMax)) {
        stop(sprintf("invalid step, must be between %d and %d", sMin, sMax))
    }
    s
}

## copy names of vector 'from' to vector or matrix 'x'
copyNames <- function(x, from) {
    if(is.null(dim(x))) {
        names(x) <- names(from)
    } else rownames(x) <- names(from)
    x
}

## copy column names of matrix 'from' to vector 'x'
copyColnames <- function(x, from) {
    names(x) <- colnames(from)
    x
}

## utility function to get default labels for plot
defaultLabels <- function(x) UseMethod("defaultLabels")
defaultLabels.seqModel <- function(x) {
    as.character(seq_along(removeIntercept(coef(x))))
}
defaultLabels.grplars <- function(x) {
    assign <- x$assign
    labels <- split(as.character(assign), assign)
    p <- sapply(labels, length)  # number of variables in each block
    append <- which(p > 1)
    labels[append] <- mapply(function(l, p) paste(l, seq_len(p), sep="."), 
        labels[append], p[append], SIMPLIFY=FALSE)
    unsplit(labels, assign)
}

## utility function to get default main plot title
defaultMain <- function() "Coefficient path" 

## find argument names of functions
#findArgNames <- function(..., removeDots = TRUE) {
#    argNames <- unique(unlist(lapply(list(...), function(f) names(formals(f)))))
#    if(removeDots) {
#        argNames <- setdiff(argNames, "...")
#    }
#    argNames
#}

## find indices of h smallest observations
findSmallest <- function(x, h) {
    # call C++ function
    callBackend <- getBackend()
    callBackend("R_findSmallest", R_x=as.numeric(x), R_h=as.integer(h)) + 1
}

## get a call function
# this returns a function that either
# 1) simply evaluates a supplied function for the basic arguments if there are
#    no additional arguments in list format
# 2) evaluates a supplied function with 'do.call' if there are additional 
#    arguments in list format
getCallFun <- function(args) {
    if(length(args) == 0) function(..., fun, args) fun(...)
    else function(..., fun, args) do.call(fun, c(list(...), args))
}

## get a component (coefficients, fitted values, residuals, ...) for certain 
## steps of the model sequence
## the component is assumed to be a matrix
getComponent <- function(x, component, s, ...) UseMethod("getComponent")
getComponent.seqModel <- function(x, component, s, ...) {
    comp <- x[[component]]      # extract the specified component
    if(missing(s)) s <- x$sOpt  # use the optimal step size as default
    if(!is.null(s)) {
        s <- checkSteps(s, sMin=1, sMax=ncol(comp))  # check steps
        comp <- comp[, s]  # extract selected steps
    }
    comp
}
getComponent.rlars <- getComponent.grplars <- function(x, component, s, ...) {
    comp <- x[[component]]      # extract the specified component
    if(missing(s)) s <- x$sOpt  # use the optimal step size as default
    if(!is.null(s)) {
        s <- checkSteps(s, sMin=0, sMax=ncol(comp)-1)  # check steps
        comp <- comp[, s + 1]  # extract selected steps
    }
    comp
}

## get the control object for model functions
#' @import robustbase MASS
getRegControl <- function(fun) {
	if(identical(fun, lmrob)) {
		fun <- .lmrob.fit
		useFormula <- FALSE
	} else if(identical(fun, rlm)) {
		fun <- .rlm
		useFormula <- FALSE
	} else useFormula <- TRUE
	list(fun=fun, useFormula=useFormula)
}

.lmrob.fit <- function(x, y, control, max.it = 500, k.max = 2500, ...) {
	if(missing(control)) {
		control <- lmrob.control(max.it=max.it, k.max=k.max, ...)
	}
	lmrob.fit(x, y, control=control)
}

.rlm <- function(x, y, maxit = 500, ...) rlm(x, y, maxit=maxit, ...)

## get steps of model sequence
getSteps <- function(x) UseMethod("getSteps")

getSteps.seqModel <- function(x) seq_along(x$df)

getSteps.rlars <- getSteps.grplars <- function(x) seq_along(x$df) - 1

## compute coefficients of hyperplane through data points
hyperplane <- function(x) {
    p <- ncol(x)
    y <- -x[, p]  # right hand side
    x <- x[, -p, drop=FALSE]
    tx <- t(x)
    theta <- solve(tx %*% x) %*% tx %*% y
    c(theta, 1)
}

## obtain the degrees of freedom of a model (number of nonzero parameters)
modelDf <- function(beta, tol = .Machine$double.eps^0.5) {
    length(which(abs(beta) > tol))
}

# find indices of h smallest observations
partialOrder <- function(x, h) {
    # call C++ function
    callBackend <- getBackend()
    callBackend("R_partialOrder", R_x=as.numeric(x), R_h=as.integer(h)) + 1
}

## find indices of h smallest observations
#partialSort <- function(x, h) {
#    # call C++ function
#    callBackend <- getBackend()
#    callBackend("R_partialSort", R_x=as.numeric(x), R_h=as.integer(h))
#}

## prepend something to column names of a matrix
prependcolnames <- function(x, prefix) {
    colnames(x) <- paste(prefix, colnames(x), sep=".")
    x
}
