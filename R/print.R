# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method print seqModel
print.seqModel <- function(x, zeros = TRUE, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print coefficients of optimal submodel
    cat("\nCoefficients of optimal submodel:\n")
    print(coef(x, zeros=zeros), ...)
    # print index of optimal submodel
    cat(sprintf("\nIndex of optimal submodel: %d\n", x$sOpt))
    # return object invisibly
    invisible(x)
}

#' @S3method print grplars
print.grplars <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print LARS sequence of groups
    active <- x$active
    steps <- seq_along(active)
    active <- t(active)
    dimnames(active) <- list("Group", steps)
    cat("\nSequence of LARS moves:\n")
    print(active, ...)
    # print coefficients of optimal LARS submodel
    cat("\nCoefficients of optimal LARS submodel:\n")
    print(coef(x, zeros=FALSE), ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print tslarsP
print.tslarsP <- function(x, ...) {
    # print "grplars" model
    print.grplars(x, ...)
    # print lag length
#    cat("\nLag length:\n")
#    print(x$p, ...)
    cat(sprintf("\nLag length: %d\n", x$p))
    # return object invisibly
    invisible(x)
}

#' @S3method print tslars
print.tslars <- function(x, ...) {
    # print "grplars" model with optimal lag length
    pOpt <- x$pOpt
    xOpt <- x$pFit[[pOpt]]
    xOpt$call <- x$call
    print.grplars(xOpt, ...)
    # print optimal lag length
#    cat("\nOptimal lag length:\n")
#    print(pOpt, ...)
    cat(sprintf("\nOptimal lag length: %d\n", pOpt))
    # return object invisibly
    invisible(x)
}

#' @S3method print rlars
print.rlars <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print LARS sequence of groups
    active <- x$active
    steps <- seq_along(active)
    active <- t(active)
    dimnames(active) <- list("Var", steps)
    cat("\nSequence of LARS moves:\n")
    print(active, ...)
    # print coefficients of optimal LARS submodel
    cat("\nCoefficients of optimal LARS submodel:\n")
    print(coef(x, zeros=FALSE), ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print sparseLTS
print.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"), 
        zeros = TRUE, ...) {
    # initializations
    fit <- match.arg(fit)
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print coefficients of sparse LTS model
    cat("\nCoefficients:\n")
    coefficients <- coef(x, fit=fit, zeros=zeros)
    print(coefficients, ...)
    # print penalty parameter and robust scale estimate
    text <- c("Penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
        lambda <- as.matrix(x$lambda)
        dimnames(lambda) <- list(text[1], "")
        print(lambda, ...)
        scale <- t(c(x$scale, x$raw.scale))
        dimnames(scale) <- list(text[2], colnames(coefficients))
        print(scale, ...)
    } else {
        scale <- if(fit == "reweighted") x$scale else x$raw.scale
        info <- matrix(c(x$lambda, scale), 2, 1)
        dimnames(info) <- list(text, "")
        print(info, ...)
    }
    # return object invisibly
    invisible(x)
}

#' @S3method print sparseLTSGrid
print.sparseLTSGrid <- function(x, fit = c("reweighted", "raw", "both"), 
        zeros = TRUE, ...) {
    # initializations
    fit <- match.arg(fit)
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print coefficients of sparse LTS model
    cat("\nCoefficients of optimal submodel:\n")
    coefficients <- coef(x, fit=fit, zeros=zeros)
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        colnames(coefficients) <- fits
    }
    print(coefficients, ...)
    # print penalty parameter and robust scale estimate
    text <- c("Optimal penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
        sOpt <- c(x$sOpt, x$raw.sOpt)
        lambda <- x$lambda[sOpt]
        scale <- c(x$scale[sOpt[1]], x$raw.scale[sOpt[2]])
        info <- rbind(lambda, scale)
        dimnames(info) <- list(text, colnames(coefficients))
        cat("\n")
    } else {
        if(fit == "reweighted") {
            sOpt <- x$sOpt
            scale <- x$scale
        } else {
            sOpt <- x$raw.sOpt
            scale <- x$raw.scale
        }
        info <- matrix(c(x$lambda[sOpt], scale[sOpt]), 2, 1)
        dimnames(info) <- list(text, "")
    }
    print(info, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print cvSeqModel
#' @import cvTools
print.cvSeqModel <- function(x, ...) {
#    # print function call
#    if(!is.null(call <- x$call)) {
#        cat("\nCall:\n")
#        dput(x$call)
#    }
    # print cross-validation results
    cvTools:::print.cvSelect(x, best=FALSE, ...)
    # print optimal model if requested
    bestFit <- x$cv[x$best, "Fit"]
    if(is.factor(bestFit)) bestFit <- as.character(bestFit)
    bestFit <- as.numeric(bestFit)
    if(length(bestFit) > 1) {
        names(bestFit) <- names(x$best)
        cat("\nOptimal step:\n")
    } else {
        bestFit <- matrix(bestFit, dimnames=list("Optimal step:", ""))
    }
    print(bestFit, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print cvSparseLTS
#' @import cvTools
print.cvSparseLTS <- function(x, ...) {
#    # print function call
#    if(!is.null(call <- x$call)) {
#        cat("\nCall:\n")
#        dput(x$call)
#    }
    # print cross-validation results
    cvTools:::print.cvTuning(x, best=FALSE, ...)
    # print optimal value for penalty parameter
    optimalLambda <- x$tuning[x$best, "lambda"]
    if(length(optimalLambda) > 1) {
        names(optimalLambda) <- names(x$best)
        cat("\nOptimal lambda:\n")
    } else {
        optimalLambda <- matrix(optimalLambda, 
            dimnames=list("Optimal lambda:", ""))
    }
    print(optimalLambda, ...)
    # return object invisibly
    invisible(x)
}
