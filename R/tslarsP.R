# ----------------------------------------------
# Author: Andreas Alfons
#         KU Leuven
#
# based on code by: Sarah Gelper
#                   Erasmus University Rotterdam
# ----------------------------------------------

#' Least angle regression for time series data with fixed lag length
#' 
#' Sequence groups of candidate predictors and their respective lagged values 
#' according to their predictive content and find the optimal model along the 
#' sequence.  Note that lagged values of the response are included as a 
#' predictor group as well.
#' 
#' @aliases print.tslarsP
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{tslarsP} is called.
#' @param x  a numeric matrix or data frame containing the candidate predictor 
#' series.
#' @param y  a numeric vector containing the response series.
#' @param h  an integer giving the forecast horizon (defaults to 1).
#' @param p  an integer giving the number of lags in the model (defaults to 2).
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
#' If \code{fit} is \code{FALSE}, an integer vector containing the indices of 
#' the sequenced predictor groups.
#'  
#' Otherwise an object of class \code{"tslarsP"} (inheriting from class 
#' \code{"grplars"}) with the following components:
#' @returnItem active  an integer vector containing the sequence of predictor 
#' groups.
#' @returnItem coefficients  a numeric matrix in which each column contains the 
#' regression coefficients of the corresponding submodel along the sequence.
#' @returnItem fitted.values  a numeric matrix in which each column contains the 
#' fitted values of the corresponding submodel along the sequence.
#' @returnItem residuals  a numeric matrix in which each column contains the 
#' residuals of the corresponding submodel along the sequence.
#' @returnItem crit  a character string specifying the optimality criterion used 
#' for selecting the final model.
#' @returnItem critValues  a numeric vector containing the values of the 
#' optimality criterion from the submodels along the sequence.
#' @returnItem df  an integer vector containing the degrees of freedom of the 
#' submodels along the sequence (i.e., the number of estimated coefficients).
#' @returnItem sOpt  an integer giving the optimal submodel.
#' @returnItem muY  numeric; the mean of the response.
#' @returnItem sigmaY  numeric; the standard deviation of the response.
#' @returnItem muX  a numeric vector containing the means of the predictor 
#' variables.
#' @returnItem sigmaX  a numeric vector containing the standard deviations of 
#' the predictor variables.
#' @returnItem x  the matrix of candidate predictor series (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response series (if \code{model} is \code{TRUE}).
#' @returnItem assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @returnItem robust  a logical indicating whether a robust fit was computed 
#' (\code{FALSE} for function \code{tslarsP}).
#' @returnItem h  the forecast horizon.
#' @returnItem p  the number of lags in the model.
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
#' @seealso \code{\link{coef.seqModel}}, \code{\link{fitted.seqModel}}, 
#' \code{\link{residuals.seqModel}}, \code{\link{predict.tslarsP}}, 
#' \code{\link{plot.seqModel}}, \code{\link{tslars}}, \code{\link{rtslarsP}}
#' 
#' @keywords regression ts
#' 
#' @export

tslarsP <- function(x, ...) UseMethod("tslarsP")


#' @rdname tslarsP
#' @method tslarsP formula
#' @export

tslarsP.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("tslarsP")
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
    out <- tslarsP.default(x, y, ...)
    if(inherits(out, "grplars")) {
        out$call <- call  # add call to return object
        out$terms <- mt   # add model terms to return object
    }
    out
}


#' @rdname tslarsP
#' @method tslarsP default
#' @export

tslarsP.default <- function(x, y, h = 1, p = 2, sMax = NA, fit = TRUE, 
        crit = "BIC", model = TRUE, ...) {
    ## call fit function with classical functions for center, scale, 
    ## correlation and regression
    call <- match.call()  # get function call
    call[[1]] <- as.name("tslarsP")
    out <- tslarsPFit(x, y, h=h, p=p, sMax=sMax, robust=FALSE, centerFun=mean, 
        scaleFun=sd, fit=fit, crit=crit, model=model)
    if(inherits(out, "grplars")) out$call <- call  # add call to return object
    out
}


#' Robust least angle regression for time series data with fixed lag length
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
#' \code{rtslarsP} is called.
#' @param x  a numeric matrix or data frame containing the candidate predictor 
#' series.
#' @param y  a numeric vector containing the response series.
#' @param h  an integer giving the forecast horizon (defaults to 1).
#' @param p  an integer giving the number of lags in the model (defaults to 2).
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
#' If \code{fit} is \code{FALSE}, an integer vector containing the indices of 
#' the sequenced predictor groups.
#'  
#' Otherwise an object of class \code{"tslarsP"} (inheriting from class 
#' \code{"grplars"}) with the following components:
#' @returnItem active  an integer vector containing the sequence of predictor 
#' groups.
#' @returnItem coefficients  a numeric matrix in which each column contains the 
#' regression coefficients of the corresponding submodel along the sequence.
#' @returnItem fitted.values  a numeric matrix in which each column contains the 
#' fitted values of the corresponding submodel along the sequence.
#' @returnItem residuals  a numeric matrix in which each column contains the 
#' residuals of the corresponding submodel along the sequence.
#' @returnItem crit  a character string specifying the optimality criterion used 
#' for selecting the final model.
#' @returnItem critValues  a numeric vector containing the values of the 
#' optimality criterion from the submodels along the sequence.
#' @returnItem df  an integer vector containing the degrees of freedom of the 
#' submodels along the sequence (i.e., the number of estimated coefficients).
#' @returnItem sOpt  an integer giving the optimal submodel.
#' @returnItem muY  numeric; the mean of the cleaned response.
#' @returnItem sigmaY  numeric; the standard deviation of the cleaned response.
#' @returnItem muX  a numeric vector containing the means of the cleaned 
#' predictor variables.
#' @returnItem sigmaX  a numeric vector containing the standard deviations of 
#' the cleaned predictor variables.
#' @returnItem x  the matrix of candidate predictor series (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response series (if \code{model} is \code{TRUE}).
#' @returnItem assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @returnItem robust  a logical indicating whether a robust fit was computed 
#' (\code{TRUE} for function \code{rtslarsP}).
#' @returnItem w  a numeric vector giving the data cleaning weights.
#' @returnItem h  the forecast horizon.
#' @returnItem p  the number of lags in the model.
#' @returnItem call  the matched function call.
#' 
#' @note The predictor group of lagged values of the response is indicated by 
#' the index 0.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{coef.seqModel}}, \code{\link{fitted.seqModel}}, 
#' \code{\link{residuals.seqModel}}, \code{\link{predict.tslarsP}}, 
#' \code{\link{plot.seqModel}}, \code{\link{rtslars}}, \code{\link{tslarsP}}
#' 
#' @keywords regression robust ts
#' 
#' @export

rtslarsP <- function(x, ...) UseMethod("rtslarsP")


#' @rdname rtslarsP
#' @method rtslarsP formula
#' @export

rtslarsP.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rtslarsP")
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
    out <- rtslarsP.default(x, y, ...)
    if(inherits(out, "grplars")) {
        out$call <- call  # add call to return object
        out$terms <- mt   # add model terms to return object
    }
    out
}


#' @rdname rtslarsP
#' @method rtslarsP default
#' @export

rtslarsP.default <- function(x, y, h = 1, p = 2, sMax = NA, 
        centerFun = median, scaleFun = mad, regFun = lmrob, 
        regArgs = list(), winsorize = FALSE, const = 2, prob = 0.95, 
        combine = c("min", "euclidean", "mahalanobis"), ncores = 1, 
        cl = NULL, fit = TRUE, crit = "BIC", model = TRUE, ...) {
    ## call fit function with classical functions for center, scale, 
    ## correlation and regression
    call <- match.call()  # get function call
    call[[1]] <- as.name("rtslarsP")
    out <- tslarsPFit(x, y, h=h, p=p, sMax=sMax, robust=TRUE, 
        centerFun=centerFun, scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, 
        winsorize=winsorize, const=const, prob=prob, combine=combine, 
        ncores=ncores, cl=cl, fit=fit, crit=crit, model=model)
    if(inherits(out, "grplars")) out$call <- call  # add call to return object
    out
}


## fit function for fixed lag length that allows to specify functions for 
## center, scale, correlation and regression
tslarsPFit <- function(x, y, 
    ## arguments for time series
    h = 1,      # forecast horizon 
    p = 2,      # maximum number of lags
    sMax = NA,  # number of predictors to be ranked
    ## arguments for correlation and short regression estimates
    robust = FALSE,     # logical indicating whether methods are robust
    centerFun = mean,   # center function
    scaleFun = sd,      # scale function
    regFun = lm.fit,    # (short) regression function
    regArgs = list(),   # additional arguments for (short) regression function
    winsorize = FALSE,  # logical indicating whether data should be winsorized
    const = 2,          # tuning constant for winsorization
    prob = 0.95,        # tuning constant for winsorization
    combine = c("min", "euclidean", "mahalanobis"),
    ncores = 1,    # number of cores for parallel computing
    cl = NULL,     # cluster for parallel computing
    ## arguments for optimal model selection
    fit = TRUE,    # logical indicating whether to fit models along sequence
    crit = "BIC",  # character string specifying the optimality criterion
    model = TRUE
) {
    ## initializations
    n <- length(y)
    x <- as.matrix(x)
    if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
    m <- ncol(x)
    assign <- rep(seq_len(m+1), each=p)
	model <- isTRUE(model)
    ## call workhorse function and modify return object
    out <- grplarsWork(fitBlocks(x, y, h, p), y[(p+h):n], sMax=sMax, 
		assign=assign, dummy=FALSE, robust=robust, centerFun=centerFun, 
		scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, winsorize=winsorize, 
		const=const, prob=prob, combine=combine, ncores=ncores, cl=cl, fit=fit, 
        crit=crit, model=model)
    # modify return object
    if(inherits(out, "grplars")) {
        out$active <- out$active - 1  # lagged response should have index 0
        if(model) {
            # include original data rather than derived data
            out$x <- x
            out$y <- y
        }
        out$assign <- out$assign - 1  # lagged response should have index 0
        out$h <- h
        out$p <- p
        class(out) <- c("tslarsP", class(out))  # "tslarsP" inherits from "grplars"
    } else out <- out - 1  # lagged response should have index 0
    out
}
