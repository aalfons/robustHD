# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## workhorse function for groupwise LARS
#' @import Rcpp 
#' @import RcppArmadillo
#' @import parallel

grplarsInternal <- function(x, y, sMax = NA, assign, dummy = TRUE, 
        robust = FALSE, centerFun = mean, scaleFun = sd, regFun = lm.fit, 
        regArgs = list(), winsorize = FALSE, const = 2, prob = 0.95, 
        combine = c("min", "euclidean", "mahalanobis"), fit = TRUE, 
        crit = c("BIC", "PE"), splits = foldControl(), cost = rmspe, 
        costArgs = list(), selectBest = c("hastie", "min"), seFactor = 1, 
        ncores = 1, cl = NULL, seed = NULL, model = TRUE) {
    ## initializations
    n <- length(y)
    assignList <- split(seq_len(length(assign)), assign)  # column indices for each block in list form
    m <- length(assignList)  # number of blocks
    p <- sapply(assignList, length)  # number of variables in each block
    adjust <- length(unique(p)) > 1  # adjust for block length?
    robust <- isTRUE(robust)
    if(robust) {
        dummy <- sapply(dummy, isTRUE)
        haveDummies <- any(dummy)
        regControl <- getRegControl(regFun)
        regFun <- regControl$fun  # if possible, do not use formula interface
        callRegFun <- getCallFun(regArgs)
        winsorize <- isTRUE(winsorize) && !haveDummies  # do not winsorize if there are dummies
        combine <- match.arg(combine)
    }
    fit <- isTRUE(fit)
    if(!is.null(seed)) set.seed(seed)
    if(is.na(ncores)) ncores <- detectCores()  # use all available cores
    if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
        ncores <- 1  # use default value
        warning("invalid value of 'ncores'; using default value")
    } else ncores <- as.integer(ncores)
    if(fit || (robust && !winsorize)) {
        # check whether parallel computing should be used
        haveCl <- inherits(cl, "cluster")
        haveNcores <- !haveCl && ncores > 1
        useParallel <- haveNcores || haveCl
        # set up multicore or snow cluster if not supplied
        if(haveNcores) {
            if(.Platform$OS.type == "windows") {
                cl <- makePSOCKcluster(rep.int("localhost", ncores))
            } else cl <- makeForkCluster(ncores)
            on.exit(stopCluster(cl))
        }
        if(useParallel) {
            # set seed of the random number stream
            if(!is.null(seed)) clusterSetRNGStream(cl, iseed=seed)
            else if(haveNcores) clusterSetRNGStream(cl)
        }
    }
    addModel <- isTRUE(model)
    if(robust) {
        if(haveDummies) {
            # use standardization with mean/SD for dummies 
            # and with centerFun/scaleFun otherwise
            z <- robStandardize(y, centerFun, scaleFun)  # robustly standardize response
            dummy <- rep(dummy, length.out=sum(p))
            xs <- robStandardizeDummy(x, dummy, centerFun, scaleFun)
        } else {
            z <- robStandardize(y, centerFun, scaleFun)  # robustly standardize response
            xs <- robStandardize(x, centerFun, scaleFun)
        }
        muY <- attr(z, "center")
        sigmaY <- attr(z, "scale")
        muX <- attr(xs, "center")
        sigmaX <- attr(xs, "scale")
        if(winsorize) {
            if(is.null(const)) const <- 2
            w <- winsorize(cbind(y, x), standardized=TRUE, 
                const=const, prob=prob, return="weights")
        } else {
            # clean data in a limited sense: there may still be correlation 
            # outliers between the blocks, but these should not be a problem
            # compute weights from robust regression for each block
            if(combine == "min") {
                # define function to compute weights for each predictor group
                getWeights <- function(i, x, y) {
                    x <- x[, i, drop=FALSE]
                    if(regControl$useFormula) {
                        fit <- callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
                    } else fit <- callRegFun(x, y, fun=regFun, args=regArgs)
                    sqrt(weights(fit))
                }
                # compute weights for each predictor group
                if(useParallel) {
                    w <- parSapply(cl, assignList, getWeights, xs, z)
                } else w <- sapply(assignList, getWeights, xs, z)
                w <- apply(w, 1, min)  # take smallest weight for each observation
                # observations can have zero weight, in which case the number 
                # of observations needs to be adjusted
                n <- length(which(w > 0))
            } else {
                # define function to compute scaled residuals for each 
                # predictor group
                getResiduals <- function(i, x, y) {
                    x <- x[, i, drop=FALSE]
                    if(regControl$useFormula) {
                        fit <- callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
                    } else fit <- callRegFun(x, y, fun=regFun, args=regArgs)
                    residuals <- residuals(fit)
                    sigma <- fit$scale
                    if(is.null(sigma)) sigma <- fit$s
                    if(is.null(sigma)) sigma <- scaleFun(residuals)
                    residuals/sigma
                }
                # compute scaled residuals for each predictor group
                if(useParallel) {
                    residuals <- parSapply(cl, assignList, getResiduals, xs, z)
                } else residuals <- sapply(assignList, getResiduals, xs, z)
                # obtain weights from scaled residuals
                if(combine == "euclidean") {
                    # assume diagonal structure of the residual correlation matrix
                    # and compute weights based on resulting mahalanobis distances
                    d <- qchisq(prob, df=m)  # quantile of the chi-squared distribution
                    w <- pmin(sqrt(d/rowSums(residuals^2)), 1)
                } else {
                    # get weights from multivariate winsorization of residuals
                    w <- winsorize(residuals, standardized=TRUE, 
                        const=const, prob=prob, return="weights")
                }
            }
        }
        z <- standardize(w*z)  # standardize cleaned response
        xs <- standardize(w*xs)
        # center and scale of response
        muY <- muY + attr(z, "center")
        sigmaY <- sigmaY * attr(z, "scale")
        # center and scale of candidate predictor variables
        muX <- muX + attr(xs, "center")
        sigmaX <- sigmaX * attr(xs, "scale")
    } else {
        z <- standardize(y)   # standardize response
        xs <- standardize(x)
        # center and scale of response
        muY <- attr(z, "center")
        sigmaY <- attr(z, "scale")
        # center and scale of candidate predictor variables
        muX <- attr(xs, "center")
        sigmaX <- attr(xs, "scale")
    }
    sMax <- checkSMax(sMax, n, m)  # check maximum number of steps
    
    ## call C++ function
    active <- .Call("R_fastGrplars", R_x=xs, R_y=z, R_sMax=as.integer(sMax[1]), 
        R_assign=assignList, R_ncores=ncores, PACKAGE="robustHD") + 1
    
    ## choose optimal model according to specified criterion
    if(fit) {
        # add ones to matrix of predictors to account for intercept
        x <- addIntercept(x)
        # call function to fit models along the sequence
        s <- if(is.na(sMax[2])) NULL else 0:sMax[2]
        out <- seqModel(x, y, active=active, s=s, assign=assignList, 
            robust=robust, regFun=regFun, useFormula=regControl$useFormula, 
            regArgs=regArgs, crit=crit, splits=splits, cost=cost, 
            costArgs=costArgs, selectBest=selectBest, seFactor=seFactor, 
            cl=cl)
        # add center and scale estimates
        out$muY <- muY
        out$sigmaY <- sigmaY
        out$muX <- muX
        out$sigmaX <- sigmaX
        if(addModel) {
            # add model data to result
            out$x <- x
            out$y <- y
        }
        out$assign <- assign  # this is helpful for some methods
        out$robust <- robust
        if(robust) {
            out$w <- w
#            out$centerFun <- centerFun
#            out$scaleFun <- scaleFun
#            out$winsorize <- winsorize
#            out$const <- const
#            out$prob <- prob
#            out$combine <- combine
        }
        class(out) <- c("grplars", class(out))
        out
    } else active
}


## find possible step sizes for groupwise LARS by solving quadratic equation
#' @export
findStepSizes <- function(r, a, corY, corU, tau) {
    mapply(function(corY, corU, tau, r, a) {
            # quadratic equation to be solved
            comp <- c(r^2 - corY^2, 2 * (corU*corY - a*r), a^2 - tau^2)
            # solution of quadratic equation
            gamma <- Re(polyroot(comp))
            min(gamma[gamma >= 0])
        }, corY, corU, tau, MoreArgs=list(r=r, a=a))
}
