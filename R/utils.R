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
    if(!is.numeric(s) || length(s) == 0 || any(!is.finite(s)) || 
            any(s < sMin) || any(s > sMax)) {
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

defaultLabels.sparseLTSGrid <- defaultLabels.seqModel

## utility function to get default main plot title
defaultMain <- function() "Coefficient path" 

## drop dimension in case of matrix with one column
dropCol <- function(x) {
    d <- dim(x)
    if(is.null(d[2]) || d[2] != 1) x
    else if(d[1] == 1) {
        # drop() drops all names for a 1x1 matrix
        names <- rownames(x)
        x <- drop(x)
        names(x) <- names
        x
    } else drop(x)
}

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

## construct blocks of original and lagged values for time series models
fitBlocks <- function(x, y, h = 1, p = 2, intercept = FALSE) {
    n <- length(y)
	tsBlocks(x, y, p=p, subset=-((n-h+1):n), intercept=intercept)
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

getSteps.seqModel <- function(x) seq_along(x$df) - 1

getSteps.sparseLTSGrid <- function(x) seq_along(x$lambda)

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

## construct blocks of original and lagged values for prediction from time 
## series models
newdataBlocks <- function(x, y, h = 1, p = 2, intercept = TRUE) {
    n <- length(y)
    tsBlocks(x, y, p=p, subset=(n-h-p+2):n, intercept=intercept)
}

## find indices of h smallest observations
partialOrder <- function(x, h) {
    # call C++ function
    callBackend <- getBackend()
    callBackend("R_partialOrder", R_x=as.numeric(x), R_h=as.integer(h)) + 1
}

## find indices of h smallest observations
#partialOrder <- function(x, h) {
#    # call C++ function
#    .CallBackend <- getBackendEnv("backend")
#    .CallBackend("R_partialOrder", R_x=as.numeric(x), R_h=as.integer(h)) + 1
#}

## find indices of h smallest observations
#partialSort <- function(x, h, backend) {
#    # call C++ function
#    .CallBackend <- getBackendEnv("backend")
#    .CallBackend("R_partialSort", R_x=as.numeric(x), R_h=as.integer(h))
#}

## get pca scores corresponding to eigenvalues larger than 1
pcaScores <- function(x, kMax) {
    # check maximum number of principal components
    d <- dim(x)
    kMax <- rep(kMax, length.out=1)
    if(!isTRUE(is.finite(kMax)) || !isTRUE(kMax <= min(d[2], d[1]-1))) {
        kMax <- min(d[2], d[1]-1)
    }
    # fit PCA and extract scores
    pca <- PCAgrid(x, k=kMax, scores=TRUE)
    sdev <- pca$sdev
    k <- which.min(sdev[sdev >= 1])
    pca$scores[, seq_len(k), drop=FALSE]
}

## prepend something to column names of a matrix
prependColnames <- function(x, prefix) {
    colnames(x) <- paste(prefix, colnames(x), sep=".")
    x
}
