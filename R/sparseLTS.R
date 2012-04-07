# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Sparse least trimmed squares regression
#' 
#' Compute least trimmed squares regression with an \eqn{L_{1}}{L1} penalty on 
#' the regression coefficients, which allows for sparse model estimates.
#' 
#' @aliases print.sparseLTS
#' 
#' @param formula  a formula describing the model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{sparseLTS} is called.
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a non-negative numeric value giving the penalty parameter.
#' @param mode  a character string specifying the type of penalty parameter.  If 
#' \code{"lambda"}, \code{lambda} gives penalty parameter directly.  If 
#' \code{"fraction"}, the smallest value of the penalty parameter that sets all 
#' coefficients to 0 is first estimated based on bivariate winsorization, then 
#' \code{lambda} gives the fraction of that estimate to be used (hence 
#' \code{lambda} should be in the interval [0,1] in that case).
#' @param alpha  a numeric value giving the percentage of the residuals for 
#' which the \eqn{L_{1}}{L1} penalized sum of squares should be minimized (the 
#' default is 0.75).
#' @param intercept  a logical indicating whether a constant term should be 
#' included in the model (the default is \code{TRUE}).
#' @param nsamp  a numeric vector giving the number of subsamples to be used in 
#' the two phases of the algorithm.  The first element gives the number of 
#' initial subsamples to be used.  The second element gives the number of 
#' subsamples to keep after the first phase of \code{ncstep} C-steps.  For 
#' those remaining subsets, additional C-steps are performed until 
#' convergence.  The default is to first perform \code{ncstep} C-steps on 500 
#' initial subsamples, and then to keep the 10 subsamples with the lowest value 
#' of the objective function for additional C-steps until convergence.
#' @param initial  a character string specifying the type of initial subsamples 
#' to be used.  If \code{"sparse"}, the lasso fit given by three randomly 
#' selected data points is first computed.  The corresponding initial subsample 
#' is then formed by the fraction \code{alpha} of data points with the smallest 
#' squared residuals.  Note that this is optimal from a robustness point of 
#' view, as the probability of including an outlier in the initial lasso fit is 
#' minimized.  If \code{"hyperplane"}, a hyperplane through \eqn{p} randomly 
#' selected data points is first computed, where \eqn{p} denotes the number of 
#' variables.  The corresponding initial subsample is then again formed by the 
#' fraction \code{alpha} of data points with the smallest squared residuals.  
#' Note that this cannot be applied if \eqn{p} is larger than the number of 
#' observations.  Nevertheless, the probability of including an outlier 
#' increases with increasing dimension \eqn{p}.  If \code{"random"}, the 
#' initial subsamples are given by a fraction \code{alpha} of randomly 
#' selected data points.  Note that this leads to the largest probability of 
#' including an outlier.
#' @param ncstep  a positive integer giving the number of C-steps to perform on 
#' all subsamples in the first phase of the algorithm (the default is to 
#' perform two C-steps).
#' @param use.correction  currently ignored.  Small sample correction factors 
#' may be added in the future.
#' @param tol  a small positive numeric value giving the tolerance for 
#' convergence.
#' @param eps  a small positive numeric value used to determine whether the 
#' variability within a variable is too small (an effective zero).
#' @param use.Gram  a logical indicating whether the Gram matrix of the 
#' explanatory variables should be precomputed in the lasso fits (the default 
#' is \code{TRUE}).  If the number of variables is large (e.g., larger than the 
#' number of observations), computation may be faster when this is set to 
#' \code{FALSE}.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param model  a logical indicating whether the data \code{x} and \code{y} 
#' should be added to the return object.  If \code{intercept} is \code{TRUE}, 
#' a column of ones is added to \code{x} to account for the intercept.
#' @param \dots  additional arguments to be passed down.
#' 
#' @returnClass sparseLTS
#' @returnItem best  an integer vector giving the best subset of \eqn{h} 
#' observations found and used for computing the raw estimates.
#' @returnItem objective  numeric; the value of the sparse LTS objective 
#' function, i.e., the \eqn{L_{1}}{L1} penalized sum of the \eqn{h} smallest 
#' squared residuals from the raw fit.
#' @returnItem coefficients  a numeric vector of coefficient estimates of the 
#' reweighted fit (including the intercept if \code{intercept} is \code{TRUE}).
#' @returnItem fitted.values  a numeric vector containing the fitted values of 
#' the response of the reweighted fit.
#' @returnItem residuals  a numeric vector containing the residuals of the 
#' reweighted fit.
#' @returnItem center  a numeric value giving the robust center estimate of the 
#' reweighted residuals.
#' @returnItem scale  a numeric value giving the robust scale estimate of the 
#' reweighted residuals.
#' @returnItem lambda  a numeric value giving the penalty parameter.
#' @returnItem intercept  a logical indicating whether the model includes a 
#' constant term.
#' @returnItem alpha  a numeric value giving the percentage of the residuals for 
#' which the \eqn{L_{1}}{L1} penalized sum of squares was minimized.
#' @returnItem quan  the number \eqn{h} of observations used to compute the raw 
#' estimates.
#' @returnItem cnp2  a numeric value giving the consistency factor applied to 
#' the scale estimate of the reweighted residuals.
#' @returnItem weights  an integer vector containing binary weights that 
#' indicate outliers, i.e., the weights are \eqn{1} for observations with 
#' reasonably small reweighted residuals and \eqn{0} for observations with 
#' large reweighted residuals.
#' @returnItem df  an integer giving the degrees of freedom of the obtained 
#' reweighted model fit, i.e., the number of nonzero coefficient estimates.
#' @returnItem raw.coefficients  a numeric vector of coefficient estimates of 
#' the raw fit (including the intercept if \code{intercept} is \code{TRUE}).  
#' @returnItem raw.residuals  a numeric vector containing the residuals of 
#' the raw fit.
#' @returnItem raw.center  a numeric value giving the robust center estimate of 
#' the raw residuals.
#' @returnItem raw.scale  a numeric value giving the robust scale estimate of 
#' the raw residuals.
#' @returnItem raw.cnp2  a numeric value giving the consistency factor applied 
#' to the scale estimate of the raw residuals.
#' @returnItem raw.weights  an integer vector containing binary weights that 
#' indicate outliers of the raw fit, i.e., the weights used for the reweighted 
#' fit.
#' @returnItem x  the predictor matrix (if \code{model} is \code{TRUE}).
#' @returnItem y  the response variable (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{sparseLTSGrid}}, \code{\link{coef.sparseLTS}}, 
#' \code{\link{fitted.sparseLTS}}, \code{\link{plot.sparseLTS}}, 
#' \code{\link{predict.sparseLTS}}, \code{\link{residuals.sparseLTS}}, 
#' \code{\link{weights.sparseLTS}}, \code{\link[robustbase]{ltsReg}}
#' 
#' @example inst/doc/examples/example-sparseLTS.R
#' 
#' @keywords regression robust
#' 
#' @export 
#' @import Rcpp 
#' @import RcppEigen
#' @useDynLib robustHD

sparseLTS <- function(x, ...) UseMethod("sparseLTS")


#' @rdname sparseLTS
#' @method sparseLTS formula
#' @export

sparseLTS.formula <- function(formula, data, ...) {
    # get function call
    call <- match.call()
    call[[1]] <- as.name("sparseLTS")
    # prepare model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if(is.empty.model(mt)) stop("empty model")
    # extract response and candidate predictors from model frame
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    # check if the specified model contains an intercept
    # if so, remove the column for intercept and use 'intercept=TRUE'
    # otherwise use 'intercept=FALSE'
    whichIntercept <- match("(Intercept)", colnames(x), nomatch = 0)
    intercept <- whichIntercept > 0
    if(intercept) x <- x[, -whichIntercept, drop = FALSE]
    # call default method
    fit <- sparseLTS(x, y, intercept=intercept, ...)
    fit$call <- call  # add call to return object
    fit$terms <- mt   # add model terms to return object
    fit
}


#' @rdname sparseLTS
#' @method sparseLTS default
#' @export

sparseLTS.default <- function(x, y, lambda, mode = c("lambda", "fraction"), 
        alpha = 0.75, intercept = TRUE, nsamp = c(500, 10), 
        initial = c("sparse", "hyperplane", "random"), ncstep = 2, 
        use.correction = TRUE, tol = .Machine$double.eps^0.5, 
        eps = .Machine$double.eps, use.Gram = TRUE, seed = NULL, 
        model = TRUE, ...) {
    ## initializations
    call <- match.call()
    call[[1]] <- as.name("sparseLTS")
    n <- length(y)
    x <- addColnames(as.matrix(x))
    d <- dim(x)
    if(!isTRUE(n == d[1])) stop(sprintf("'x' must have %d rows", n))
    if(missing(lambda)) {
        # if penalty parameter is not supplied, use a small fraction of a 
        # robust estimate of the smallest value that sets all coefficients 
        # to zero
        lambda <- 0.05
        mode <- "fraction"
    } else {
        # otherwise check the supplied penalty parameter
        lambda <- rep(lambda, length.out=1)
        if(!is.numeric(lambda) || !is.finite(lambda)) {
            stop("missing or invalid value of 'lambda'")
        }
        if(lambda < 0) {
            lambda <- 0
            warning("negative value for 'lambda', using no penalization")
        }
        mode <- match.arg(mode)
    }
    intercept <- isTRUE(intercept)
    if(mode == "fraction" && lambda > 0) { 
        # fraction of a robust estimate of the smallest value for the penalty 
        # parameter that sets all coefficients to zero (based on bivariate 
        # winsorization)
        lambda <- lambda * 
            lambda0(x, y, intercept=intercept, tol=tol, eps=eps, ...)
    }
    alpha <- rep(alpha, length.out=1)
    if(!isTRUE(is.numeric(alpha) && 0.5 <= alpha && alpha <= 1)) {
        stop("'alpha' must be between 0.5 and 1")
    }
    h <- ceiling(alpha*n)  # subset sizes are different from 'ltsReg'
    nsamp <- rep(nsamp, length.out=2)
    if(!is.numeric(nsamp) || any(!is.finite(nsamp))) {
        nsamp <- c(500, 10)
        warning("missing or infinite values in 'nsamp'; using default values")
    }
    ncstep <- rep(ncstep, length.out=1)
    if(!is.numeric(ncstep) || !is.finite(ncstep)) {
        ncstep <- 2
        warning("missing or infinite value of 'ncstep'; using default value")
    }
    tol <- rep(tol, length.out=1)
    if(!is.numeric(tol) || !is.finite(tol)) {
        tol <- .Machine$double.eps^0.5
        warning("missing or infinite value of 'tol'; using default value")
    }
    eps <- rep(eps, length.out=1)
    if(!is.numeric(eps) || !is.finite(eps)) {
        eps <- .Machine$double.eps
        warning("missing or infinite value of 'eps'; using default value")
    }
    use.Gram <- isTRUE(use.Gram)
    if(!is.null(seed)) set.seed(seed)
    
    ## compute initial subsets
#   subsets <- randomSubsets(nrow(x), h, nsamp=nsamp[1]) - 1
    subsets <- initialSubsets(x, y, h, nsamp[1], type=initial, lambda=lambda, 
        intercept=intercept, eps=eps, use.Gram=use.Gram) - 1
    
    ## call C++ function
    fit <- .Call("R_fastSparseLTS", R_x=x, R_y=y, R_lambda=lambda, 
        R_subsets=subsets, R_intercept=intercept, 
        R_ncstep=as.integer(ncstep), R_nkeep=as.integer(nsamp[2]), 
        R_tol=tol, R_eps=eps, R_useGram=use.Gram, PACKAGE="robustHD")
    best <- sort.int(fit$best + 1)
    objective <- fit$crit
    
    ## correct scale estimate if requested
    q <- qnorm(0.9875)  # quantile of the normal distribution
    qn <- qnorm((h+n)/ (2*n))  # quantile for consistency factor
    cdelta <- 1 / sqrt(1 - (2*n)/(h/qn) * dnorm(qn))  # consistency factor
    s <- fit$scale * cdelta
    ## compute 0/1 weights identifying outliers
    ok <- abs((fit$residuals - fit$center)/s) <= q  # good observations
    raw.weights <- as.integer(ok)
    
    ## compute reweighted estimator
    # keep information on raw estimator
    raw.fit <- fit
    raw.cdelta <- cdelta
    raw.s <- s
    nOk <- sum(raw.weights)  # number of good observations
    # compute reweighted estimate
    fit <- fastLasso(x, y, lambda=lambda, subset=which(ok), 
        intercept=intercept, eps=eps, use.Gram=use.Gram)
#    # -----
#    # experiment: use only the selected variables for the lasso in the 
#    #             reweighting step
#    raw.coef <- raw.fit$coefficients
#    select <- which(abs(raw.coef[-1]) > tol)  # selected variables
#    fit <- fastLasso(x[, select, drop=FALSE], y, lambda=lambda, 
#        subset=which(ok), intercept=intercept, eps=eps, use.Gram=use.Gram)
#    coef <- rep.int(0, length(raw.coef))
#    if(intercept) select <- c(1, select + 1)
#    coef[select] <- fit$coefficients
#    # modify object
#    fit$coefficients <- coef
#    # -----
    if(intercept) x <- addIntercept(x)  # add column to account for intercept
    if(nOk < n) {
        qn <- qnorm((nOk+n)/ (2*n))  # quantile for consistency factor
        cdelta <- 1 / sqrt(1-(2*n)/(nOk/qn)*dnorm(qn))  # consistency factor
    } else cdelta <- 1  # consistency factor not necessary
    center <- sum(raw.weights*fit$residuals)/nOk
    centeredResiduals <- fit$residuals - center
    s <- sqrt(sum(raw.weights*centeredResiduals^2)/(nOk-1)) * cdelta
    weights <- as.integer(abs(centeredResiduals/s) <= q)
    
    ## compute degrees of freedom (number of nonzero parameters)
    df <- modelDf(fit$coefficients, tol)
    
    ## return object
    fit <- list(best=best, objective=objective, 
        coefficients=copyColnames(fit$coefficients, x), 
        fitted.values=copyNames(fit$fitted.values, y), 
        residuals=copyNames(fit$residuals, y), center=center, scale=s, 
        lambda=lambda, intercept=intercept, alpha=alpha, quan=h, 
        cnp2=cdelta, weights=weights, df=df, 
        raw.coefficients=copyColnames(raw.fit$coefficients, x), 
        raw.residuals=copyNames(raw.fit$residuals, y), 
        raw.center=raw.fit$center, raw.scale=raw.s, raw.cnp2=raw.cdelta, 
        raw.weights=raw.weights)
    if(isTRUE(model)) {
        fit$x <- x
        fit$y <- y
    }
    fit$call <- call
    class(fit) <- "sparseLTS"
    fit
}


#' Sparse least trimmed squares regression
#' 
#' Compute least trimmed squares regression with an \eqn{L_{1}}{L1} penalty on 
#' the regression coefficients, which allows for sparse model estimates, over 
#' a grid of values for the penalty parameter.
#' 
#' @aliases print.sparseLTSGrid
#' 
#' @param formula  a formula describing the model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{sparseLTS} is called.
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a numeric vector of non-negative numeric values to be used as 
#' penalty parameter.
#' @param mode  a character string specifying the type of penalty parameter.  If 
#' \code{"lambda"}, \code{lambda} gives the grid of values for the penalty 
#' parameter directly.  If \code{"fraction"}, the smallest value of the penalty 
#' parameter that sets all coefficients to 0 is first estimated based on 
#' bivariate winsorization, then \code{lambda} gives the fractions of that 
#' estimate to be used (hence all values of \code{lambda} should be in the 
#' interval [0,1] in that case).
#' @param crit  a character string specifying the optimality criterion to be 
#' used for selecting the final model.  Currently, only \code{"BIC"} for the 
#' Bayes information criterion is implemented.
#' @param \dots  additional arguments to be passed to \code{\link{sparseLTS}}.
#' @param model  a logical indicating whether the data \code{x} and \code{y} 
#' should be added to the return object.  If \code{intercept} is \code{TRUE}, 
#' a column of ones is added to \code{x} to account for the intercept.
#' 
#' @return An object of class \code{"sparseLTSGrid"} (inheriting from class 
#' \code{"seqModel"}) with the following components:
#' @returnItem best  an integer matrix in which each column contains the best 
#' subset of \eqn{h} observations found and used for computing the raw 
#' estimates with the corresponding penalty parameter.
#' @returnItem objective  a numeric vector giving the values of the sparse LTS 
#' objective function, i.e., the \eqn{L_{1}}{L1} penalized sum of the \eqn{h} 
#' smallest squared residuals from the raw fits.
#' @returnItem coefficients  a numeric matrix in which each column contains the 
#' coefficient estimates of the corresponding reweighted fit (including the 
#' intercept if \code{intercept} is \code{TRUE}).
#' @returnItem fitted.values  a numeric matrix in which each column contains 
#' the fitted values of the response of the corresponding reweighted fit.
#' @returnItem residuals  a numeric matrix in which each column contains 
#' the residuals of the response of the corresponding reweighted fit.
#' @returnItem center  a numeric vector giving the robust center estimates 
#' of the residuals from the reweighted fits.
#' @returnItem scale  a numeric vector giving the robust scale estimates of 
#' the residuals from the reweighted fits.
#' @returnItem lambda  a numeric vector giving the values of the penalty 
#' parameter.
#' @returnItem intercept  a logical indicating whether the model includes a 
#' constant term.
#' @returnItem alpha  a numeric value giving the percentage of the residuals for 
#' which the \eqn{L_{1}}{L1} penalized sum of squares was minimized.
#' @returnItem quan  the number \eqn{h} of observations used to compute the raw 
#' estimates.
#' @returnItem cnp2  a numeric vector giving the consistency factors applied to 
#' the scale estimates of the residuals from the reweighted fits.
#' @returnItem weights  an integer matrix in which each column contains binary 
#' weights that indicate outliers from the corresponding reweighted fit, i.e., 
#' the weights are \eqn{1} for observations with reasonably small reweighted 
#' residuals and \eqn{0} for observations with large reweighted residuals.
#' @returnItem df  an integer vector giving the degrees of freedom of the 
#' obtained reweighted model fits, i.e., the number of nonzero coefficient 
#' estimates.
#' @returnItem raw.coefficients  a numeric matrix in which each column contains 
#' the coefficient estimates of the corresponding raw fit (including the 
#' intercept if \code{intercept} is \code{TRUE}).  
#' @returnItem raw.residuals  a numeric matrix in which each column contains 
#' the residuals of the corresponding raw fit.
#' @returnItem raw.center  a numeric vector giving the robust center estimates 
#' of the residuals from the raw fits.
#' @returnItem raw.scale  a numeric vector giving the robust scale estimates of 
#' the residuals from the raw fits.
#' @returnItem raw.cnp2  a numeric vector giving the consistency factors 
#' applied to the scale estimates of the residuals from the raw fits.
#' @returnItem raw.weights  an integer matrix in which each column contains 
#' binary weights that indicate outliers of the corresponding raw fit, i.e., 
#' the weights used for the reweighted fits.
#' @returnItem crit  a character string specifying the optimality criterion used 
#' for selecting the optimal model.
#' @returnItem critValues  a numeric vector containing the values of the 
#' optimality criterion from the reweighted fits.
#' @returnItem sOpt  an integer giving the optimal reweighted fit.
#' @returnItem raw.critValues  a numeric vector containing the values of the 
#' optimality criterion from the raw fits.
#' @returnItem raw.sOpt  an integer giving the optimal raw fit.
#' @returnItem x  the predictor matrix (if \code{model} is \code{TRUE}).
#' @returnItem y  the response variable (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{sparseLTS}}, 
#' \code{\link[=coef.sparseLTS]{coef.sparseLTSGrid}}, 
#' \code{\link[=fitted.sparseLTS]{fitted.sparseLTSGrid}}, 
#' \code{\link[=plot.seqModel]{plot.sparseLTSGrid}}, 
#' \code{\link{diagnosticPlot}}, 
#' \code{\link[=predict.sparseLTS]{predict.sparseLTSGrid}}, 
#' \code{\link[=residuals.sparseLTS]{residuals.sparseLTSGrid}}, 
#' \code{\link[=weights.sparseLTS]{weights.sparseLTSGrid}}, 
#' 
#' @example inst/doc/examples/example-sparseLTSGrid.R
#' 
#' @keywords multivariate robust
#' 
#' @export 
#' @import Rcpp 
#' @import RcppEigen
#' @useDynLib robustHD

sparseLTSGrid <- function(x, ...) UseMethod("sparseLTSGrid")


#' @rdname sparseLTSGrid
#' @method sparseLTSGrid formula
#' @export

sparseLTSGrid.formula <- function(formula, data, ...) {
    # get function call
    call <- match.call()
    call[[1]] <- as.name("sparseLTSGrid")
    # prepare model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if(is.empty.model(mt)) stop("empty model")
    # extract response and candidate predictors from model frame
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    # check if the specified model contains an intercept
    # if so, remove the column for intercept and use 'intercept=TRUE'
    # otherwise use 'intercept=FALSE'
    whichIntercept <- match("(Intercept)", colnames(x), nomatch = 0)
    intercept <- whichIntercept > 0
    if(intercept) x <- x[, -whichIntercept, drop = FALSE]
    # call default method
    fit <- sparseLTSGrid(x, y, intercept=intercept, ...)
    fit$call <- call  # add call to return object
    fit$terms <- mt   # add model terms to return object
    fit
}


#' @rdname sparseLTSGrid
#' @method sparseLTSGrid default
#' @export

sparseLTSGrid.default <- function(x, y, lambda, mode = c("lambda", "fraction"), 
        crit = "BIC", ..., model = TRUE) {
    ## initializations
    call <- match.call()
    call[[1]] <- as.name("sparseLTSGrid")
    if(missing(lambda)) {
        # if penalty parameter is not supplied, use a small fraction of a 
        # robust estimate of the smallest value that sets all coefficients 
        # to zero
        lower <- if(nrow(x) > ncol(x)) 0 else 0.1
        lambda <- seq(from=1, to=lower, by=-0.1)
        mode <- "fraction"
    } else {
        # otherwise check the supplied penalty parameter
        if(!is.numeric(lambda) || length(lambda) == 0 || any(!is.finite(lambda))) {
            stop("missing or invalid value of 'lambda'")
        }
        if(any(negative <- lambda < 0)) {
            lambda[negative] <- 0
            warning("negative value for 'lambda', using no penalization")
        }
        mode <- match.arg(mode)
        lambda <- sort(unique(lambda), decreasing=TRUE)
    }
    if(mode == "fraction" && any(lambda > 0)) { 
        # fraction of a robust estimate of the smallest value for the penalty 
        # parameter that sets all coefficients to zero (based on bivariate 
        # winsorization)
        lambda <- lambda * lambda0(x, y, ...)
    }
    crit <- match.arg(crit)
    # fit sparse LTS models along supplied grid
    fit <- lapply(lambda, 
        function(l, ...) {
            sparseLTS(x, y, lambda=l, mode="lambda", ..., model=FALSE)
        }, ...)
    # select the optimal reweighted and raw model via BIC
    if(crit == "BIC") {
        critValues <- sapply(fit, BIC, fit="both")
        raw.critValues <- critValues["raw",]
        critValues <- critValues["reweighted",]
        sOpt <- which.min(critValues)
        raw.sOpt <- which.min(raw.critValues)
    }
    # combine information from the models into suitable data structures
    names(fit) <- seq_along(lambda)
    best <- sapply(fit, function(x) x$best)
    objective <- sapply(fit, function(x) x$objective)
    coef <- sapply(fit, coef, fit="reweighted")
    fitted <- sapply(fit, fitted, fit="reweighted")
    residuals <- sapply(fit, residuals, fit="reweighted")
    center <- sapply(fit, function(x) x$center)
    scale <- sapply(fit, function(x) x$scale)
    lambda <- sapply(fit, function(x) x$lambda)
    intercept <- fit[[1]]$intercept
    alpha <- fit[[1]]$alpha
    quan <- fit[[1]]$quan
    cnp2 <- sapply(fit, function(x) x$cnp2)
    weights <- sapply(fit, weights, fit="reweighted")
    df <- sapply(fit, function(x) x$df)
    raw.coef <- sapply(fit, coef, fit="raw")
    raw.residuals <- sapply(fit, residuals, fit="raw")
    raw.center <- sapply(fit, function(x) x$raw.center)
    raw.scale <- sapply(fit, function(x) x$raw.scale)
    raw.cnp2 <- sapply(fit, function(x) x$raw.cnp2)
    raw.weights <- sapply(fit, weights, fit="raw")
    # construct return object
    fit <- list(best=best, objective=objective, coefficients=coef, 
        fitted.values=fitted, residuals=residuals, center=center, 
        scale=scale, lambda=lambda, intercept=intercept, alpha=alpha, 
        quan=quan, cnp2=cnp2, weights=weights, df=df, 
        raw.coefficients=raw.coef, raw.residuals=raw.residuals, 
        raw.center=raw.center, raw.scale=raw.scale, raw.cnp2=raw.cnp2, 
        raw.weights=raw.weights, crit=crit, critValues=critValues, 
        sOpt=sOpt, raw.critValues=raw.critValues, raw.sOpt=raw.sOpt)
    if(isTRUE(model)) {
        if(intercept) x <- addIntercept(x)
        fit$x <- x
        fit$y <- y
    }
    fit$call <- call
    class(fit) <- c("sparseLTSGrid", "seqModel")
    fit
}
