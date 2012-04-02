# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

fitModels <- function(
    x,  # matrix of predictors
    y,  # response
    s = NULL,  # integer vector giving the steps of the submodels to be fitted
    assign = NULL,       # integer vector giving the group assignment of the variables
    robust = FALSE,      # logical indicating whether methods are robust
    regFun = lm.fit,     # regression function
    useFormula = FALSE,  # logical indicating whether to use formula interface
    regArgs = list(),    # additional arguments for regression function
    active,              # integer vector giving the sequence of predictors
    final = TRUE,        # logical indicating whether to return final object
    crit = "BIC",        # character string specifying the optimality criterion
    class                # character string giving the class of the return object
) {
    n <- length(y)
    haveAssign <- !is.null(assign)
    if(haveAssign && !is.list(assign)) {
        # column indices for each block in list form
        assign <- split(seq_len(length(assign)), assign)
    }
    if(robust) callRegFun <- getCallFun(regArgs)
#    # add ones to matrix of predictors to account for intercept
#    x <- addIntercept(x, check=TRUE)
    if(haveAssign) {
        # the default for 's' is to fit models as long as there are twice as 
        # many observations as predictors
        # if 's' is supplied, keep only those steps for which there are more 
        # observations than predictors
        if(is.null(s)) {
            sMax <- floor(n/2)
            s <- 0:min(length(active), sMax)
            dfMax <- sMax + 1
        } else dfMax <- n
        # compute degrees of freedom of the submodels along sequence
        # R uses (df + 1) in BIC penalty for (weighted) linear models!!!
        whichMax <- which.max(s)
        firstActive <- head(active, s[whichMax])
        p <- sapply(assign[firstActive], length)  # number of variables in each block
        df <- unname(cumsum(c(1, p)))[s+1]        # degrees of freedom
        # only fit submodels while the degrees of freedom does not become 
        # larger than the requested maximum
        if(df[whichMax] > dfMax) {
            keep <- which(df <= dfMax)
            s <- s[keep]
            whichMax <- which.max(s)
            firstActive <- head(active, s[whichMax])
            df <- df[keep]
        }
        # blockwise sequenced variables (including intercept)
        sequenced <- c(1, unlist(assign[firstActive], use.names=FALSE) + 1)
    } else {
        # the default for 's' is to fit models as long as there are twice as 
        # many observations as predictors
        if(is.null(s)) s <- 0:min(length(active), floor(n/2))
        # compute degrees of freedom of the submodels along sequence
        # R uses (df + 1) in BIC penalty for (weighted) linear models!!!
        df <- s + 1  # account for intercept
        # sequenced variables (including intercept)
        sequenced <- c(1, head(active, max(s)) + 1)
    }
    # fit submodels
    # number of variables to use is one less than degrees of freedom
    if(robust) {
        if(useFormula) {
            models <- lapply(df, 
                function(k) {
                    x <- x[, sequenced[seq_len(k)], drop=FALSE]
                    callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
                })
        } else {
            models <- lapply(df, 
                function(k) {
                    x <- x[, sequenced[seq_len(k)], drop=FALSE]
                    callRegFun(x, y, fun=regFun, args=regArgs)
                })
        }
    } else {
        models <- lapply(df, 
            function(k) lm.fit(x[, sequenced[seq_len(k)], drop=FALSE], y))
    }
    # construct matrix of coefficents
    coef <- matrix(0, nrow=ncol(x), ncol=length(s), 
        dimnames=list(colnames(x), s))
    for(k in seq_len(ncol(coef))) {
        coef[sequenced[seq_len(df[k])], k] <- coef(models[[k]])
    }
    if(final) {
        # compute fitted values and residuals along sequence
        fitted <- sapply(models, fitted)
        residuals <- sapply(models, residuals)
        colnames(fitted) <- colnames(residuals) <- s
        # compute optimality criterion along sequence
        if(crit == "BIC") {
            # R uses (df + 1) in BIC penalty for (weighted) linear models!!!
            penalty <- log(n)  # penalty for BIC
            if(robust) {
                # BIC based on residual scale estimate from robust regression
                critFun <- function(k) AIC(models[[k]], k=penalty)
            } else {
                critFun <- function(k) {
                    n * (log(2 * pi) + 1 - log(n) + 
                            log(sum(residuals[, k]^2))) + df[k] * penalty
                }
            }
        }
        critValues <- sapply(seq_len(ncol(residuals)), critFun)
        if(crit == "BIC") sOpt <- which.min(critValues)
        sOpt <- sOpt - 1
        ## construct return object
        out <- list(active=active, coefficients=coef, fitted.values=fitted, 
            residuals=residuals, crit=crit, critValues=critValues, df=df, 
            sOpt=sOpt)
    } else {
        ## construct return object
        # this is only used for cross-validation, hence only coefficents 
        # are included in the object
        out <- list(coefficients=coef)
    }
    class(out) <- c(class, "seqModel")
    out
}
