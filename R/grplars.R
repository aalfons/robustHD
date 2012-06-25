# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Groupwise least angle regression
#' 
#' Sequence groups of candidate predictors according to their predictive 
#' content and find the optimal model along the sequence.
#' 
#' @aliases print.grplars findStepSizes
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{grplars} is called.
#' @param x  a matrix or data frame containing the candidate predictors.
#' @param y  a numeric vector containing the response.
#' @param sMax  an integer vector of length two.  If a single integer is 
#' supplied, it is recycled.  The first element gives the number of predictor 
#' groups to be sequenced.  If it is \code{NA} (the default), predictor 
#' groups are sequenced as long as there are no singularity issues.  The second 
#' element gives the maximum number of predictor groups to be included in the 
#' final model.  If it is \code{NA} (the default), predictor groups may be 
#' added to the model as long as there are twice as many observations as 
#' predictor variables.
#' @param assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
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
#' Otherwise an object of class \code{"grplars"} (inheriting from class 
#' \code{"seqModel"}) with the following components:
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
#' @returnItem x  the matrix of candidate predictors (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response (if \code{model} is \code{TRUE}).
#' @returnItem assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @returnItem robust  a logical indicating whether a robust fit was computed 
#' (\code{FALSE} for function \code{grplars}).
#' @returnItem call  the matched function call.
#' 
#' @note \code{findStepSizes} is a utility function that computes the step size 
#' for each inactive predictor group.  It is only exported so it can be called 
#' by the underlying C++ code for sequencing the predictor groups.  Hence it is 
#' not expected to be called by the user and not documented.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{coef.seqModel}}, \code{\link{fitted.seqModel}}, 
#' \code{\link{residuals.seqModel}}, \code{\link{predict.seqModel}}, 
#' \code{\link{plot.seqModel}}, \code{\link{rgrplars}}
#' 
#' @keywords regression
#' 
#' @export

grplars <- function(x, ...) UseMethod("grplars")


#' @rdname grplars
#' @method grplars formula
#' @export

grplars.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("grplars")
    # prepare model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 1  # ensure model with intercept
#    if(is.empty.model(mt)) stop("empty model")
    # extract response and candidate predictors from model frame
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    ## call default method
    out <- grplars.default(x, y, ...)
    if(inherits(out, "grplars")) {
        out$call <- call  # add call to return object
        out$terms <- mt   # add model terms to return object
    }
    out
}


#' @rdname grplars
#' @method grplars data.frame
#' @export

grplars.data.frame <- function(x, y, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("grplars")
    x <- model.matrix(~ ., data=x)   # convert data.frame to design matrix
    ## call default method
    out <- grplars.default(x, y, ...)
    if(inherits(out, "grplars")) out$call <- call  # add call to return object
    out
}


#' @rdname grplars
#' @method grplars default
#' @export

grplars.default <- function(x, y, sMax = NA, assign, fit = TRUE, 
        crit = "BIC", ncores = 1, model = TRUE, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("grplars")
    # if argument 'assign' is not supplied, check if predictor matrix has 
    # attribute "assign" (as generated by model.matrix())
    if(missing(assign)) {
        assign <- attr(x, "assign")
        if(is.null(assign)) {
            # take each predictor to be its own group
            assign <- seq_len(ncol(x))
        }
    }
    # check if the predictor matrix contains column for intercept and 
    # remove it if necessary
    if(isTRUE(assign[1] == 0)) {
        x <- removeIntercept(x)
        assign <- assign[-1]
    }
    ## call fit function with classical functions for center, scale, 
    ## correlation and regression
    out <- grplarsFit(x, y, sMax=sMax, assign=assign, robust=FALSE, 
        centerFun=mean, scaleFun=sd, fit=fit, crit=crit, ncores=ncores, 
        model=model)
    if(inherits(out, "grplars")) out$call <- call  # add call to return object
    out
}



#' Robust groupwise least angle regression
#' 
#' Robustly sequence groups of candidate predictors according to their 
#' predictive content and find the optimal model along the sequence.
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{rgrplars} is called.
#' @param x  a matrix or data frame containing the candidate predictors.
#' @param y  a numeric vector containing the response.
#' @param sMax  an integer vector of length two.  If a single integer is 
#' supplied, it is recycled.  The first element gives the number of predictor 
#' groups to be sequenced.  If it is \code{NA} (the default), predictor 
#' groups are sequenced as long as there are no singularity issues.  The second 
#' element gives the maximum number of predictor groups to be included in the 
#' final model.  If it is \code{NA} (the default), predictor groups may be 
#' added to the model as long as there are twice as many observations as 
#' predictor variables.
#' @param assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @param dummy  a logical vector indicating whether the predictors are dummy 
#' variables.
#' @param centerFun  a function to compute a robust estimate for the center 
#' (defaults to \code{\link[stats]{median}}).
#' @param scaleFun  a function to compute a robust estimate for the scale 
#' (defaults to \code{\link[stats]{mad}}).
#' @param regFun  a function to compute robust linear regressions that can be 
#' interpreted as weighted least squares (defaults to 
#' \code{\link[robustbase]{lmrob}}).
#' @param regArgs  a list of arguments to be passed to \code{regFun}.
#' @param combine  a character string specifying how to combine the data 
#' cleaning weights from the robust regressions with each predictor group.  
#' Possible values are \code{"min"} for taking the minimum weight for each 
#' observation, \code{"euclidean"} for weights based on Euclidean distances 
#' of the multivariate set of standardized residuals (i.e., multivariate 
#' winsorization of the standardized residuals assuming independence), or 
#' \code{"mahalanobis"} for weights based on Mahalanobis distances of the 
#' multivariate set of standardized residuals (i.e., multivariate winsorization 
#' of the standardized residuals).
#' @param const  numeric; tuning constant for multivariate winsorization to be 
#' used in the initial corralation estimates based on adjusted univariate 
#' winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in multivariate 
#' winsorization (defaults to 0.95).
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
#' Otherwise an object of class \code{"grplars"} (inheriting from class 
#' \code{"seqModel"}) with the following components:
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
#' @returnItem x  the matrix of candidate predictors (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response (if \code{model} is \code{TRUE}).
#' @returnItem assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @returnItem robust  a logical indicating whether a robust fit was computed 
#' (\code{TRUE} for function \code{rgrplars}).
#' @returnItem w  a numeric vector giving the data cleaning weights.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{coef.seqModel}}, \code{\link{fitted.seqModel}}, 
#' \code{\link{residuals.seqModel}}, \code{\link{predict.seqModel}}, 
#' \code{\link{plot.seqModel}}, \code{\link{grplars}}
#' 
#' @keywords regression
#' 
#' @export

rgrplars <- function(x, ...) UseMethod("rgrplars")


#' @rdname rgrplars
#' @method rgrplars formula
#' @export

rgrplars.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rgrplars")
    # prepare model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 1  # ensure model with intercept
#    if(is.empty.model(mt)) stop("empty model")
    # extract response and candidate predictors from model frame
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    ## call wrapper around default method
    out <- rgrplarsDefault(x, y, mt, ...)
    if(inherits(out, "grplars")) {
        out$call <- call  # add call to return object
        out$terms <- mt   # add model terms to return object
    }
    out
}


#' @rdname rgrplars
#' @method rgrplars data.frame
#' @export

rgrplars.data.frame <- function(x, y, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rgrplars")
    # convert data.frame to design matrix
    mf <- model.frame(~ ., data=x)
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    ## call default method
    out <- rgrplarsDefault(x, y, mt, ...)
    if(inherits(out, "grplars")) out$call <- call  # add call to return object
    out
}


#' @rdname rgrplars
#' @method rgrplars default
#' @export

rgrplars.default <- function(x, y, sMax = NA, assign, dummy, 
        centerFun = median, scaleFun = mad, regFun = lmrob, regArgs = list(), 
        combine = c("min", "euclidean", "mahalanobis"), const = 2, prob = 0.95, 
        ncores = 1, cl = NULL, fit = TRUE, crit = "BIC", model = TRUE, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rgrplars")
    # if argument 'assign' is not supplied, check if predictor matrix has 
    # attribute "assign" (as generated by model.matrix())
    if(missing(assign)) {
        assign <- attr(x, "assign")
        if(is.null(assign)) {
            # take each predictor to be its own group
            assign <- seq_len(ncol(x))
        }
    }
    # check if the predictor matrix contains column for intercept and 
    # remove it if necessary
    if(isTRUE(assign[1] == 0)) {
        x <- removeIntercept(x)
        assign <- assign[-1]
    }
    # if argument 'dummy' is not supplied, check which columns in the predictor 
    # matrix are are dummies
    if(missing(dummy)) dummy <- apply(x, 2, function(x) all(x %in% c(0, 1)))
    ## call fit function with supplied functions for center, scale, 
    ## correlation and regression
    out <- grplarsFit(x, y, sMax=sMax, assign=assign, dummy=dummy, robust=TRUE, 
        centerFun=centerFun, scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, 
		combine=combine, const=const, prob=prob, ncores=ncores, cl=cl, fit=fit, 
        crit=crit, model=model)
    if(inherits(out, "grplars")) out$call <- call  # add call to return object
    out
}

rgrplarsDefault <- function(x, y, mt, assign, dummy, ...) {
    haveIntercept <- as.logical(attr(mt, "intercept"))
    # extract group assignment of the variables
    if(missing(assign)) {
        assign <- attr(x, "assign")
        if(haveIntercept) assign <- assign[-1]
    }
    # check which columns in the predictor matrix are are dummies
    if(missing(dummy)) {
        dummy <- attr(mt, "term.labels") %in% names(attr(x, "contrasts"))
        dummy <- dummy[assign]
    }
    ## call default method
    if(haveIntercept) x <- removeIntercept(x) # remove column for intercept
    rgrplars.default(x, y, assign=assign, dummy=dummy, ...)
}


## fit function that allows to specify functions for center, scale, correlation 
## and regression
grplarsFit <- function(x, y, 
    ## general arguments
    sMax = NA,     # number of predictors to be ranked
    assign,        # integer vector giving the group assignment of the variables 
    dummy = TRUE,  # logical vector indicating dummy variables
    ## arguments for scaling, correlation and short regression estimates
    robust = FALSE,     # logical indicating whether methods are robust
    centerFun = mean,   # center function
    scaleFun = sd,      # scale function
	regFun = lm.fit,    # (short) regression function
    regArgs = list(),   # additional arguments for (short) regression function
	combine = c("min", "euclidean", "mahalanobis"),
	const = 2,     # tuning constant for winsorization
	prob = 0.95,   # tuning constant for winsorization
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
    x <- addColnames(as.matrix(x))
    if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
    ## call workhorse function
    grplarsWork(x, y, sMax=sMax, assign=assign, dummy=dummy, robust=robust, 
        centerFun=centerFun, scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, 
        const=const, prob=prob, combine=combine, ncores=ncores, cl=cl, fit=fit, 
        crit=crit, model=model)
}
