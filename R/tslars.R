# ----------------------------------------------
# Author: Andreas Alfons
#         KU Leuven
#
# based on code by: Sarah Gelper
#                   Erasmus University Rotterdam
# ----------------------------------------------

#' Least angle regression for time series data
#' 
#' Sequence groups of candidate predictors and their respective lagged values 
#' according to their predictive content and find the optimal model along the 
#' sequence.  Note that lagged values of the response are included as a 
#' predictor group as well.
#' 
#' @aliases print.tslars
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{tslars} is called.
#' @param x  a numeric matrix or data frame containing the candidate predictor 
#' series.
#' @param y  a numeric vector containing the response series.
#' @param h  an integer giving the forecast horizon (defaults to 1).
#' @param pMax  an integer giving the maximum number of lags in the model 
#' (defaults to 3).
#' @param sMax  an integer vector of length two.  If a single integer is 
#' supplied, it is recycled.  The first element gives the number of predictor 
#' groups to be sequenced.  If it is \code{NA} (the default), predictor 
#' groups are sequenced as long as there are no singularity issues.  The second 
#' element gives the maximum number of predictor groups to be included in the 
#' final model.  If it is \code{NA} (the default), predictor groups may be 
#' added to the model as long as there are twice as many observations as 
#' predictor variables.
#' @param fit  a logical indicating whether to fit submodels along the sequence 
#' (\code{TRUE}, the default) or to simply return the sequence (\code{FALSE}).
#' @param crit  a character string specifying the optimality criterion to be 
#' used for selecting the final model.  Currently, only \code{"BIC"} for the 
#' Bayes information criterion is implemented.
#' @param model  a logical indicating whether the model data should be included 
#' in the returned object.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return 
#' If \code{fit} is \code{FALSE}, an integer matrix in which each column 
#' contains the indices of the sequenced predictor groups for the corresponding 
#' lag length.
#'  
#' Otherwise an object of class \code{"tslars"} with the following components:
#' @returnItem pFit  a list containing the fits for the respective lag lengths 
#' (see \code{\link{tslarsP}}).
#' @returnItem pOpt  an integer giving the optimal number of lags.
#' @returnItem pMax  the maximum number of lags considered.
#' @returnItem x  the matrix of candidate predictor series (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response series (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#' 
#' @note The predictor group of lagged values of the response is indicated by 
#' the index 0.
#' 
#' @author Andreas Alfons, based on code by Sarah Gelper
#' 
#' @references 
#' Gelper, S. and Croux, C. (2010) Time series least angle regression for 
#' selecting predictive economic sentiment series. Working paper.
#' 
#' @seealso \code{\link{coef.tslars}}, \code{\link{fitted.tslars}}, 
#' \code{\link{residuals.tslars}}, \code{\link{predict.tslars}}, 
#' \code{\link{plot.tslars}}, \code{\link{tslarsP}}, \code{\link{rtslars}}
#' 
#' @keywords regression ts
#' 
#' @export

tslars <- function(x, ...) UseMethod("tslars")


#' @rdname tslars
#' @method tslars formula
#' @export

tslars.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("tslars")
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
    # remove first column for intercept, if existing
    if(attr(mt, "intercept")) x <- x[, -1, drop=FALSE]
    ## call default method
    out <- tslars.default(x, y, ...)
    if(inherits(out, "tslars")) {
        out$call <- call  # add call to return object
        out$terms <- mt   # add model terms to return object
    }
    out
}


#' @rdname tslars
#' @method tslars default
#' @export

tslars.default <- function(x, y, h = 1, pMax = 3, sMax = NA, fit = TRUE, 
        crit = "BIC", model = TRUE, ...) {
    ## call fit function with classical functions for center, scale, 
    ## correlation and regression
    call <- match.call()  # get function call
    call[[1]] <- as.name("tslars")
    out <- tslarsFit(x, y, h=h, pMax=pMax, sMax=sMax, robust=FALSE, 
        centerFun=mean, scaleFun=sd, fit=fit, crit=crit, model=model)
    if(inherits(out, "tslars")) out$call <- call  # add call to return object
    out
}


#' Robust least angle regression for time series data
#' 
#' Robustly sequence groups of candidate predictors and their respective lagged 
#' values according to their predictive content and find the optimal model 
#' along the sequence.  Note that lagged values of the response are included as 
#' a predictor group as well.
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{rtslars} is called.
#' @param x  a numeric matrix or data frame containing the candidate predictor 
#' series.
#' @param y  a numeric vector containing the response series.
#' @param h  an integer giving the forecast horizon (defaults to 1).
#' @param pMax  an integer giving the maximum number of lags in the model 
#' (defaults to 3).
#' @param sMax  an integer vector of length two.  If a single integer is 
#' supplied, it is recycled.  The first element gives the number of predictor 
#' groups to be sequenced.  If it is \code{NA} (the default), predictor 
#' groups are sequenced as long as there are no singularity issues.  The second 
#' element gives the maximum number of predictor groups to be included in the 
#' final model.  If it is \code{NA} (the default), predictor groups may be 
#' added to the model as long as there are twice as many observations as 
#' predictor variables.
#' @param centerFun  a function to compute a robust estimate for the center 
#' (defaults to \code{\link[stats]{median}}).
#' @param scaleFun  a function to compute a robust estimate for the scale 
#' (defaults to \code{\link[stats]{mad}}).
#' @param regFun  a function to compute robust linear regressions that can be 
#' interpreted as weighted least squares (defaults to 
#' \code{\link[robustbase]{lmrob}}).
#' @param regArgs  a list of arguments to be passed to \code{regFun}.
#' @param winsorize  a logical indicating whether to clean the data by 
#' multivariate winsorization.
#' @param const  numeric; tuning constant for multivariate winsorization to be 
#' used in the initial corralation estimates based on adjusted univariate 
#' winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in multivariate 
#' winsorization (defaults to 0.95).
#' @param combine  a character string specifying how to combine the data 
#' cleaning weights from the robust regressions with each predictor group.  
#' Possible values are \code{"min"} for taking the minimum weight for each 
#' observation, \code{"euclidean"} for weights based on Euclidean distances 
#' of the multivariate set of standardized residuals (i.e., multivariate 
#' winsorization of the standardized residuals assuming independence), or 
#' \code{"mahalanobis"} for weights based on Mahalanobis distances of the 
#' multivariate set of standardized residuals (i.e., multivariate winsorization 
#' of the standardized residuals).
#' @param fit  a logical indicating whether to fit submodels along the sequence 
#' (\code{TRUE}, the default) or to simply return the sequence (\code{FALSE}).
#' @param crit  a character string specifying the optimality criterion to be 
#' used for selecting the final model.  Currently, only \code{"BIC"} for the 
#' Bayes information criterion based on a robust residual scale estimate is 
#' implemented.
#' @param model  a logical indicating whether the model data should be included 
#' in the returned object.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return 
#' If \code{fit} is \code{FALSE}, an integer matrix in which each column 
#' contains the indices of the sequenced predictor groups for the corresponding 
#' lag length.
#'  
#' Otherwise an object of class \code{"tslars"} with the following components:
#' @returnItem pFit  a list containing the fits for the respective lag lengths.
#' (see \code{\link{rtslarsP}}).
#' @returnItem pOpt  an integer giving the optimal number of lags.
#' @returnItem pMax  the maximum number of lags considered.
#' @returnItem x  the matrix of candidate predictor series (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response series (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#' 
#' @note The predictor group of lagged values of the response is indicated by 
#' the index 0.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{coef.tslars}}, \code{\link{fitted.tslars}}, 
#' \code{\link{residuals.tslars}}, \code{\link{predict.tslars}}, 
#' \code{\link{plot.tslars}}, \code{\link{rtslarsP}}, \code{\link{tslars}}
#' 
#' @keywords regression robust ts
#' 
#' @export

rtslars <- function(x, ...) UseMethod("rtslars")


#' @rdname rtslars
#' @method rtslars formula
#' @export

rtslars.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rtslars")
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
    # remove first column for intercept, if existing
    if(attr(mt, "intercept")) x <- x[, -1, drop=FALSE]
    ## call default method and modify return object
    out <- rtslars.default(x, y, ...)
    if(inherits(out, "tslars")) {
        out$call <- call  # add call to return object
        out$terms <- mt   # add model terms to return object
    }
    out
}

#' @rdname rtslars
#' @method rtslars default
#' @export

rtslars.default <- function(x, y, h = 1, pMax = 3, sMax = NA, 
        centerFun = median, scaleFun = mad, regFun = lmrob, 
        regArgs = list(), winsorize = FALSE, const = 2, prob = 0.95, 
        combine = c("min", "euclidean", "mahalanobis"), 
        fit = TRUE, crit = "BIC", model = TRUE, ...) {
    ## call fit function with classical functions for center, scale, 
    ## correlation and regression
    call <- match.call()  # get function call
    call[[1]] <- as.name("rtslars")
    out <- tslarsFit(x, y, h=h, pMax=pMax, sMax=sMax, robust=TRUE, 
        centerFun=centerFun, scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, 
		winsorize=winsorize, const=const, prob=prob, combine=combine, 
		fit=fit, crit=crit, model=model)
    if(inherits(out, "tslars")) out$call <- call  # add call to return object
    out
}


## fit function that allows to specify functions for center, scale, correlation 
## and regression
tslarsFit <- function(x, y, 
    ## arguments for time series
    h = 1,      # forecast horizon 
    pMax = 3,   # maximum number of lags
    sMax = NA,  # number of predictors to be ranked
    ## arguments for scaling, correlation and short regression estimates
    robust = FALSE,     # logical indicating whether methods are robust
    centerFun = mean,   # center function
    scaleFun = sd,      # scale function
	regFun = lm.fit,    # (short) regression function
    regArgs = list(),   # additional arguments for (short) regression function
    winsorize = FALSE,  # logical indicating whether data should be winsorized
    const = 2,     # tuning constant for winsorization
	prob = 0.95,   # tuning constant for winsorization
	combine = c("min", "euclidean", "mahalanobis"),
	## arguments for optimal model selection
    fit = TRUE,    # logical indicating whether to fit models along sequence
    crit = "BIC",  # optimality criterion
	## other arguments,
    model = TRUE  # logical indicating whether model data should be added to result
) {
    ## initializations
    n <- length(y)
    x <- as.matrix(x)
    if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
    crit <- match.arg(crit)
    ## call 'tslarsPFit()' for each number of lags and choose optimal lag length
    # TODO: check 'h' and 'pMax'
    p <- seq_len(pMax)
    out <- lapply(p, 
        function(i) {
            select <- (pMax-i+1):n  # make sure to use the same observations
            tslarsPFit(x[select, , drop=FALSE], y[select], h=h, p=i, sMax=sMax, 
                robust=robust, centerFun=centerFun, scaleFun=scaleFun, 
                regFun=regFun, regArgs=regArgs, winsorize=winsorize, 
                const=const, prob=prob, combine=combine, fit=fit, crit=crit, 
                model=FALSE)
        })
    names(out) <- p
    ## find optimal lag length
    if(isTRUE(fit)) {
        if(pMax == 1) pOpt <- 1
        else {
            critValuesOpt <- sapply(out, function(x) x$critValues[x$sOpt+1])
            if(crit == "BIC") whichOpt <- which.min(critValuesOpt)
            pOpt <- p[whichOpt]
        }
        ## construct return object
#        out <- list(pFit=out, h=h, pOpt=pOpt, pMax=pMax)
        out <- list(pFit=out, pOpt=pOpt, pMax=pMax)
        if(isTRUE(model)) {
            out$x <- x
            out$y <- y
        }
        class(out) <- "tslars"
    } else out <- do.call(cbind, out)
    out
}
