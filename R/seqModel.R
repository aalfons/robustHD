# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @import parallel

seqModel <- function(x, y, active, s = NULL, assign = NULL, robust = FALSE, 
        regFun = lm.fit, useFormula = FALSE, regArgs = list(), final = TRUE, 
        crit = c("BIC", "PE"), splits = foldControl(), cost = rtmspe, 
        costArgs = list(), selectBest = c("hastie", "min"), seFactor = 1, 
        cl = NULL) {
    # initializations
    n <- length(y)
    haveAssign <- !is.null(assign)
    if(haveAssign && !is.list(assign)) {
        # column indices for each block in list form
        assign <- split(seq_len(length(assign)), assign)
    }
    if(robust) callRegFun <- getCallFun(regArgs)
    if(final) crit <- match.arg(crit)
#    # add ones to matrix of predictors to account for intercept
#    x <- addIntercept(x, check=TRUE)
    # prepare the variable sequence and the degrees of freedom of the models
    if(haveAssign) {
        # the default for 's' is to fit models as long as there are twice 
        # as many observations as predictors
        # if 's' is supplied, keep only those steps for which there are 
        # more observations than predictors
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
        # the default for 's' is to fit models as long as there are twice 
        # as many observations as predictors
        if(is.null(s)) s <- 0:min(length(active), floor(n/2))
        # compute degrees of freedom of the submodels along sequence
        # R uses (df + 1) in BIC penalty for (weighted) linear models!!!
        df <- s + 1  # account for intercept
        # sequenced variables (including intercept)
        sequenced <- c(1, head(active, max(s)) + 1)
    }
    # define function to fit the submodels along the sequence
    if(robust) {
        if(useFormula) {
            fitFun <- function(k) {
                x <- x[, sequenced[seq_len(k)], drop=FALSE]
                callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
            }
        } else {
            fitFun <- function(k) {
                x <- x[, sequenced[seq_len(k)], drop=FALSE]
                callRegFun(x, y, fun=regFun, args=regArgs)
            }
        }
    } else {
        fitFun <- function(k) lm.fit(x[, sequenced[seq_len(k)], drop=FALSE], y)
    }
    # fit submodels and find optimal one
    if(final && crit == "PE") {
        # select the optimal submodel via prediction error
        selectBest <- match.arg(selectBest)
        critValues <- perrySeqModel(x, y, active=active, s=s, assign=assign, 
            robust=robust, regFun=regFun, useFormula=useFormula, 
            regArgs=regArgs, splits=splits, cost=cost, costArgs=costArgs, 
            selectBest=selectBest, seFactor=seFactor, cl=cl)
        # fit optimal submodel
        kOpt <- critValues$best
        model <- fitFun(df[kOpt])
        # add zero coefficients for variables not in the optimal submodel
        coef <- copyColnames(rep.int(0, ncol(x)), x) 
        coef[sequenced[seq_len(df[kOpt])]] <- coef(model)
        # construct return object
        out <- list(active=active, df=df, coefficients=coef, 
            fitted.values=fitted(model), residuals=residuals(model), 
            crit=crit, critValues=critValues)
        class(out) <- "optSeqModel"
    } else {
        # fit submodels
        # number of variables to use is one less than degrees of freedom
        if(is.null(cl)) models <- lapply(df, fitFun) 
        else models <- parLapply(cl, df, fitFun)
        # construct matrix of coefficents
        coef <- matrix(0, nrow=ncol(x), ncol=length(s), 
            dimnames=list(colnames(x), s))
        for(k in seq_len(ncol(coef))) {
            coef[sequenced[seq_len(df[k])], k] <- coef(models[[k]])
        }
        if(final) {
            # extract fitted values and residuals along sequence
            fitted <- sapply(models, fitted)
            residuals <- sapply(models, residuals)
            colnames(fitted) <- colnames(residuals) <- s
            # compute optimality criterion along sequence
            if(crit == "BIC") {
                # R uses (df + 1) in BIC penalty for (weighted) linear models!!!
                penalty <- log(n)  # penalty for BIC
                if(robust) {
                    # BIC based on robust residual scale estimate
                    critFun <- function(k) AIC(models[[k]], k=penalty)
                } else {
                    critFun <- function(k) {
                        n * (log(2 * pi) + 1 - log(n) + 
                                log(sum(residuals[, k]^2))) + df[k] * penalty
                    }
                }
                critValues <- sapply(seq_len(ncol(residuals)), critFun)
                sOpt <- which.min(critValues) - 1
                # construct return object
                out <- list(active=active, df=df, coefficients=coef, 
                    fitted.values=fitted, residuals=residuals, crit=crit, 
                    critValues=critValues, sOpt=sOpt)
            }
        } else {
            # construct return object
            # this is only used for cross-validation, hence only coefficents 
            # are included in the object
            out <- list(coefficients=coef)
        }
        class(out) <- "seqModel"
    }
    out
}


## internal function for estimating the prediction error of sequential models
perrySeqModel <- function(x, y, active, s = NULL, assign = NULL, 
        robust = FALSE, regFun = lm.fit, useFormula = FALSE, 
        regArgs = list(), splits = foldControl(), cost, costArgs = list(), 
        selectBest = c("hastie", "min"), seFactor = 1, cl = NULL) {
    ## initializations
    if(!is.null(s)) {
        s <- checkSteps(s, 0, length(active))  # check steps
        s <- sort(unique(s))  # make sure steps are in the correct order
    }
    if(missing(cost)) cost <- if(robust) rtmspe else rmspe
    ## construct call to fit models in prediction error estimation
    call <- as.call(list(seqModel, active=active, s=s, assign=assign, 
            robust=robust, final=FALSE))
    if(robust) {
        call$regFun <- regFun
        call$useFormula <- useFormula
        call$regArgs <- regArgs
    }
    ## call function perryFit() to perform prediction error estimation
    out <- perryFit(call, x=x, y=y, splits=splits, cost=cost, 
        costArgs=costArgs, cl=cl)
    ## reshape, modify and return the prediction error results
    out <- perryReshape(out, selectBest=selectBest, seFactor=seFactor)
    fits(out) <- as.numeric(as.character(fits(out)))
    class(out) <- c("perrySeqModel", class(out))
    out
}
