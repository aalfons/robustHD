# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## workhorse function for groupwise LARS

#grplarsWork <- function(
#    ## basic arguments
#    x,  # matrix of predictor blocks
#    y,  # response
#    sMax = NA,  # number of predictors to be ranked and included in the model
#    assign,     # integer vector giving the group assignment of the variables 
#    dummy = TRUE,  # logical vector indicating dummy variables
#    ## arguments for scaling, correlation and short regression estimates
#    robust = FALSE,     # logical indicating whether methods are robust
#    centerFun = mean,   # center function
#    scaleFun = sd,      # scale function
#    regFun = lm.fit,    # (short) regression function
#    regArgs = list(),   # additional arguments for (short) regression function
#    winsorize = FALSE,  # logical indicating whether data should be winsorized
#    const = 2,     # tuning constant for adjusted univariate winsorization
#    prob = 0.95,   # tuning constant for multivariate winsorization
#	combine = c("min", "euclidean", "mahalanobis"),
#    ## arguments for optimal model selection
#    fit = TRUE,    # logical indicating whether to fit models along sequence
#    crit = "BIC",  # character string specifying the optimality criterion
#    ## other arguments,
#    model = TRUE  # logical indicating whether model data should be added to result
#) {
#    ## STEP 1: initializations
#    n <- length(y)
#    assignList <- split(seq_len(length(assign)), assign)  # column indices for each block in list form
#    m <- length(assignList)  # number of blocks
#    p <- sapply(assignList, length)  # number of variables in each block
#    adjust <- length(unique(p)) > 1  # adjust for block length?
#    robust <- isTRUE(robust)
#    if(robust) {
#        dummy <- sapply(dummy, isTRUE)
#        haveDummies <- any(dummy)
#        regControl <- getRegControl(regFun)
#		regFun <- regControl$fun  # if possible, do not use formula interface
#		callRegFun <- getCallFun(regArgs)
#		winsorize <- isTRUE(winsorize) && !haveDummies  # do not winsorize if there are dummies
#		combine <- match.arg(combine)
#	}
#	crit <- match.arg(crit)
#    addModel <- isTRUE(model)
#	if(robust) {
#		if(haveDummies) {
#            # use standardization with mean/MAD for dummies 
#            # and with centerFun/scaleFun otherwise
#			z <- robStandardize(y, centerFun, scaleFun)  # robustly standardize response
#            dummy <- rep(dummy, length.out=sum(p))
##            xs <- x
##            xs[, dummy] <- standardize(x[, dummy])
##            xs[, !dummy] <- robStandardize(x[, !dummy], centerFun, scaleFun)
#            xs <- robStandardizeDummy(x, dummy, centerFun, scaleFun)
#        } else {
#			z <- robStandardize(y, centerFun, scaleFun)  # robustly standardize response
#			xs <- robStandardize(x, centerFun, scaleFun)
#		}
#		muY <- attr(z, "center")
#		sigmaY <- attr(z, "scale")
#		muX <- attr(xs, "center")
#		sigmaX <- attr(xs, "scale")
#		if(winsorize) {
#            if(is.null(const)) const <- 2
#			w <- winsorize(cbind(y, x), standardized=TRUE, 
#				const=const, prob=prob, return="weights")
#		} else {
#			# clean data in a limited sense: there may still be correlation 
#			# outliers between the blocks, but these should not be a problem
#			# compute weights from robust regression for each block
#			if(combine == "min") {
#				w <- sapply(assignList, 
#					function(i, x, y) {
#						x <- x[, i, drop=FALSE]
#						if(regControl$useFormula) {
#							fit <- callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
#						} else fit <- callRegFun(x, y, fun=regFun, args=regArgs)
#						sqrt(weights(fit))
#					}, xs, z)
#                w <- apply(w, 1, min)  # take smallest weight for each observation
#                # observations can have zero weight, in which case the number 
#                # of observations needs to be adjusted
#                n <- length(which(w > 0))
#			} else {
#				# compute scaled residuals
#				residuals <- sapply(assignList, 
#					function(i, x, y) {
#						x <- x[, i, drop=FALSE]
#						if(regControl$useFormula) {
#							fit <- callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
#						} else fit <- callRegFun(x, y, fun=regFun, args=regArgs)
#						residuals <- residuals(fit)
#						sigma <- fit$scale
#						if(is.null(sigma)) sigma <- fit$s
#						if(is.null(sigma)) sigma <- scaleFun(residuals)
#						residuals/sigma
#					}, xs, z)
#				if(combine == "euclidean") {
#					# assume diagonal structure of the residual correlation matrix
#					# and compute weights based on resulting mahalanobis distances
#					d <- qchisq(prob, df=m)  # quantile of the chi-squared distribution
#					w <- pmin(sqrt(d/rowSums(residuals^2)), 1)
#				} else {
#					# get weights from multivariate winsorization of residuals
#		            w <- winsorize(residuals, standardized=TRUE, 
#						const=const, prob=prob, return="weights")
#				}
#			}
#		}
#		z <- standardize(w*z)  # standardize cleaned response
#		xs <- standardize(w*xs)
#		# center and scale of response
#		muY <- muY + attr(z, "center")
#		sigmaY <- sigmaY * attr(z, "scale")
#		# center and scale of candidate predictor variables
#		muX <- muX + attr(xs, "center")
#		sigmaX <- sigmaX * attr(xs, "scale")
#	} else {
#		z <- standardize(y)   # standardize response
#		xs <- standardize(x)
#        # center and scale of response
#        muY <- attr(z, "center")
#        sigmaY <- attr(z, "scale")
#        # center and scale of candidate predictor variables
#        muX <- attr(xs, "center")
#        sigmaX <- attr(xs, "scale")
#    }
#    sMax <- checkSMax(sMax, n, m)  # check maximum number of steps
#    
#    ## STEP 2: find first ranked block
#    tmp <- lapply(assignList, 
#        function(i, x, y) {
#            x <- x[, i, drop=FALSE]
#            model <- lm.fit(x, y)
#            Rsq <- var(fitted(model))
#            list(Rsq=Rsq, model=model)
#        }, xs, z)
#    Rsq <- unname(sapply(tmp, function(x) x$Rsq))
#    # first active block
#    if(adjust) {
#        A <- which.max(unname(Rsq / p))  # adjust for block length
#    } else A <- which.max(Rsq)
#    Ac <- seq_len(m)[-A]      # not yet sequenced blocks
#    # extract fitted values
#    tmp <- lapply(tmp, function(x) x$model)
#    zHat <- sapply(tmp, fitted)
#    
#    ## STEP 3: update active set
#    R <- diag(1, sMax[1])
#    r <- rep.int(NA, sMax[1])
#    a <- rep.int(NA, sMax[1])
#    if(adjust) adj <- sqrt(p)  # adjustment for different block lengths
#    s <- 1
#    sigma <- c(s, rep.int(NA, sMax[1]))
#    # start iterative computations
#    for(k in seq_len(sMax[1]-1)) {
#        if(k > 1) {
#            # update fitted values for new or not yet sequenced blocks
#            fixed <- A[seq_len(k-1)]
#            zHat[, -fixed] <- (zHat[, -fixed] - gamma.k * uHat) / sigma[k]
#        }
#        # standardize fitted values for k-th block
#        zHat[, A[k]] <- standardize(zHat[, A[k]])
#        # compute the equiangular vector
#        if(k == 1) {
#            r[k] <- sqrt(Rsq)[A]
#            a[k] <- 1
#            w.k <- 1
#            u.k <- zHat[, A[k]]  # equiangular vector equals first direction
#            # adjustment for unequal block size in case of dummy variables
#            if(adjust) {
#                adj.k <- adj[A[k]]
#                r[k] <- r[k] / adj.k
#                a[k] <- a[k] / adj.k
#            }
#        } else {
#            # adjustment for unequal block size in case of dummy variables is 
#            # taken care of by update formula for r[k], and are considered in 
#            # vector q.k for computation of a[k]
#            r[k] <- (r[k-1] - gamma.k * a[k-1]) / sigma[k]
#            # compute correlations between fitted values for active blocks
#            R[k, seq_len(k-1)] <- R[seq_len(k-1), k] <- 
#                apply(zHat[, A[seq_len(k-1)], drop=FALSE], 2, cor, zHat[, A[k]])
#            # other computations according to algorithm
#            invR <- solve(R[seq_len(k), seq_len(k), drop=FALSE])
#            q.k <- if(adjust) adj[A[seq_len(k)]] else rep.int(1, k)
#            a[k] <- c(t(q.k) %*% invR %*% q.k)^(-1/2)
#            w.k <- c(a[k] * invR %*% q.k)
#            u.k <- c(zHat[, A[seq_len(k)], drop=FALSE] %*% w.k)  # equiangular vector
#        }
#        # compute step size in equiangular direction
#        tmp <- lapply(assignList[Ac], 
#            function(i, x, u) {
#                x <- x[, i, drop=FALSE]
#                lm.fit(x, u)
#            }, xs, u.k)
#        uHat <- sapply(tmp, fitted)
#        if(k == 1) {
#            r.kj <- sqrt(Rsq[Ac])
#        } else {
#            # update formula takes care of adjustment for unequal block size
#            r.kj <- sqrt(r.kj^2 - 2 * gamma.k * a.kj * r.kj + 
#                    gamma.k^2 * tau.kj^2) / sigma[k]
#        }
#        a.kj <- apply(zHat[, Ac, drop=FALSE], 2, cor, u.k)
#        tau.kj <- apply(uHat, 2, sd)
#        # adjustment for unequal block size in case of dummy variables
#        if(adjust) {
#            adj.kj <- adj[Ac]
#            if(k == 1) r.kj <- r.kj / adj.kj  # update formula takes care of adjustment otherwise
#            a.kj <- a.kj / adj.kj
#            tau.kj <- tau.kj / adj.kj
#        }
#        # compute step size by solving quadratic equation
#        gamma <- mapply(function(r.kj, a.kj, tau.kj, r.k, a.k) {
#                # quadratic equation to be solved
#                comp <- c(r.k^2-r.kj^2, 2*(a.kj*r.kj-a.k*r.k), a.k^2-tau.kj^2)
#                # solution of quadratic equation
#                gamma.j <- Re(polyroot(comp))
#                min(gamma.j[gamma.j >= 0])
#            }, r.kj, a.kj, tau.kj, MoreArgs=list(r.k=r[k], a.k=a[k]))
#        whichMin <- which.min(gamma)
#        gamma.k <- gamma[whichMin]
#        # update active set and not yet sequenced blocks
#        A <- c(A, Ac[whichMin])  # update active set
#        Ac <- Ac[-whichMin]      # update not yet sequenced blocks
#        r.kj <- r.kj[-whichMin]
#        a.kj <- a.kj[-whichMin]
#        tau.kj <- tau.kj[-whichMin]
#        # update scale of response
#        sigma[k+1] <- sqrt(1 - 2*gamma.k*r[k]/a[k] + gamma.k^2)
#        s <- s * sigma[k+1]
#    }
#    
#    ## STEP 4: choose optimal model according to specified criterion
#    if(isTRUE(fit)) {
#        # add ones to matrix of predictors to account for intercept
#        x <- addIntercept(x)
#        # call function to fit models along the sequence
#        s <- if(is.na(sMax[2])) NULL else 0:sMax[2]
#        out <- fitModels(x, y, s=s, assign=assignList, robust=robust, 
#            regFun=regFun, useFormula=regControl$useFormula, regArgs=regArgs, 
#            active=A, crit=crit, class="grplars")
#        # add center and scale estimates
#        out$muY <- muY
#        out$sigmaY <- sigmaY
#        out$muX <- muX
#        out$sigmaX <- sigmaX
#        if(addModel) {
#            # add model data to result
#            out$x <- x
#            out$y <- y
#        }
#        out$assign <- assign  # this is helpful for plot() and repCV() methods
#        out$robust <- robust
#        if(robust) {
#            out$w <- w
#            # TODO: should information on how models are fit be stored?
##            out$centerFun <- centerFun
##            out$scaleFun <- scaleFun
##            out$winsorize <- winsorize
##            out$const <- const
##            out$prob <- prob
##            out$combine <- combine
#        }
#        out
#    } else A
#}

#' @import parallel
grplarsWork <- function(
    ## basic arguments
    x,  # matrix of predictor blocks
    y,  # response
    sMax = NA,  # number of predictors to be ranked and included in the model
    assign,     # integer vector giving the group assignment of the variables 
    dummy = TRUE,  # logical vector indicating dummy variables
    ## arguments for scaling, correlation and short regression estimates
    robust = FALSE,     # logical indicating whether methods are robust
    centerFun = mean,   # center function
    scaleFun = sd,      # scale function
    regFun = lm.fit,    # (short) regression function
    regArgs = list(),   # additional arguments for (short) regression function
    winsorize = FALSE,  # logical indicating whether data should be winsorized
    const = 2,     # tuning constant for adjusted univariate winsorization
    prob = 0.95,   # tuning constant for multivariate winsorization
    combine = c("min", "euclidean", "mahalanobis"),
    ncores = 1,    # number of cores for parallel computing
    cl = NULL,     # cluster for parallel computing
    ## arguments for optimal model selection
    fit = TRUE,    # logical indicating whether to fit models along sequence
    crit = "BIC",  # character string specifying the optimality criterion
    ## other arguments,
    model = TRUE  # logical indicating whether model data should be added to result
) {
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
    crit <- match.arg(crit)
    addModel <- isTRUE(model)
    if(robust) {
        if(haveDummies) {
            # use standardization with mean/MAD for dummies 
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
            # check whether parallel computing should be used
            haveCl <- !is.null(cl)
            useMC <- useSnow <- FALSE
            if(!missing(ncores)) {
                ncores <- rep(ncores, length.out=1)
                if(!is.numeric(ncores) || !is.finite(ncores) || ncores < 1) {
                    ncores <- 1  # use default value
                    warning("invalid value of 'ncores'; using default value")
                }
                if(ncores > 1) {
                    useSnow <- .Platform$OS.type == "windows"
                    useMC <- !useSnow
                }
            } else if(haveCl) useSnow <- TRUE
            # start snow cluster if not supplied (for parallel computing on 
            # Windows systems)
            if(useSnow && !haveCl) {
                cl <- makePSOCKcluster(rep.int("localhost", ncores))
                if(RNGkind()[1] == "L'Ecuyer-CMRG") clusterSetRNGStream(cl)
            }
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
                if(useMC) {
                    w <- simplify2array(mclapply(assignList, getWeights, 
                            xs, z, mc.cores=ncores))
                } else if(useSnow) {
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
                if(useMC) {
                    residuals <- simplify2array(mclapply(assignList, 
                            getResiduals, xs, z, mc.cores=ncores))
                } else if(useSnow) {
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
            # stop snow cluster if defined within this function
            if(useSnow && !haveCl) stopCluster(cl)
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
        R_assign=assignList, PACKAGE="robustHD") + 1
    
    ## choose optimal model according to specified criterion
    if(isTRUE(fit)) {
        # add ones to matrix of predictors to account for intercept
        x <- addIntercept(x)
        # call function to fit models along the sequence
        s <- if(is.na(sMax[2])) NULL else 0:sMax[2]
        out <- fitModels(x, y, s=s, assign=assignList, robust=robust, 
            regFun=regFun, useFormula=regControl$useFormula, regArgs=regArgs, 
            active=active, crit=crit, class="grplars")
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
        out$assign <- assign  # this is helpful for plot() and repCV() methods
        out$robust <- robust
        if(robust) {
            out$w <- w
            # TODO: should information on how models are fit be stored?
#            out$centerFun <- centerFun
#            out$scaleFun <- scaleFun
#            out$winsorize <- winsorize
#            out$const <- const
#            out$prob <- prob
#            out$combine <- combine
        }
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
