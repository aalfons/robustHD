# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Resampling-based prediction error for a sequential regression model
#' 
#' Estimate the prediction error of a previously fit sequential regression 
#' model such as a robust least angle regression model via (repeated) 
#' \eqn{K}-fold cross-validation, (repeated) random splitting (also known as 
#' random subsampling or Monte Carlo cross-validation), or the bootstrap.  In 
#' each iteration, the optimal model is thereby selected from the training data 
#' and used to make predictions for the test data.
#' 
#' @method perry seqModel
#' @aliases perry.rlars
#' 
#' @param object  the model fit for which to estimate the prediction error.
#' @param splits  an object of class \code{"cvFolds"} (as returned by 
#' \code{\link[perry]{cvFolds}}) or a control object of class 
#' \code{"foldControl"} (see \code{\link[perry]{foldControl}}) defining the 
#' folds of the data for (repeated) \eqn{K}-fold cross-validation, an object of 
#' class \code{"randomSplits"} (as returned by 
#' \code{\link[perry]{randomSplits}}) or a control object of class 
#' \code{"splitControl"} (see \code{\link[perry]{splitControl}}) defining 
#' random data splits, or an object of class \code{"bootSamples"} (as returned 
#' by \code{\link[perry]{bootSamples}}) or a control object of class 
#' \code{"bootControl"} (see \code{\link[perry]{bootControl}}) defining 
#' bootstrap samples.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' vectors to be passed as its first two arguments, the first corresponding to 
#' the observed values of the response and the second to the predicted values, 
#' and must return a non-negative scalar value.  The default is to use the root 
#' mean squared prediction error for non-robust models and the root trimmed 
#' mean squared prediction error for robust models (see 
#' \code{\link[perry]{cost}}).
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed to the prediction loss 
#' function \code{cost}.
#' 
#' @returnClass perry 
#' @returnItem pe  a numeric value giving the estimated prediction error.  In 
#' case of more than one replication, this gives the average value over all 
#' replications.
#' @returnItem se  a numeric value giving the estimated standard error of 
#' the prediction loss.
#' @returnItem reps  a numeric matrix with one column that contains the 
#' estimated prediction errors from all replications.  This is only returned 
#' in case of more than one replication.
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error.
#' @returnItem y  the response.
#' @returnItem yHat  a list containing the predicted values from all 
#' replications.
#' @returnItem seed  the seed of the random number generator before estimation 
#' of the prediction error.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @note Users may prefer the wrapper functions \code{\link[perry]{repCV}}, 
#' \code{\link[perry]{repRS}} and \code{\link[perry]{bootPE}}.
#' 
#' @seealso \code{\link{rlars}}, \code{\link{predict.seqModel}}, 
#' \code{\link[perry]{perry}}, \code{\link[perry]{cost}}
#' 
#' @example inst/doc/examples/example-perry.rlars.R
#' 
#' @keywords utilities
#' 
#' @export
#' @import perry

perry.seqModel <- function(object, splits = foldControl(), 
        cost, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perry")
    # retrieve data from model fit
    if(is.null(x <- object$x) || is.null(y <- object$y)) {
        if(is.null(x)) x <- try(model.matrix(object$terms), silent=TRUE)
        if(is.null(y)) y <- try(model.response(object$terms), silent=TRUE)
        if(inherits(x, "try-error") || inherits(y, "try-error")) {
            stop("model data not available")
        }
    }
    # predictor matrix is stored with column for intercept
    x <- removeIntercept(x)
    # get default cost function depending on whether model fit is robust
    robust <- if(inherits(object, "rlars")) TRUE else object$robust
    if(missing(cost)) cost <- if(robust) rtmspe else rmspe
    ## construct call to fit models in prediction error estimation
    call <- object$call
    # if the model was fitted with formula method, 'formula' and 'data' 
    # arguments are removed from call and 'x' and 'y' are used instead
    call$formula <- NULL
    call$data <- NULL
    ## call function perryFit() to perform prediction error estimation
    out <- perryFit(call, x=x, y=y, splits=splits, cost=cost, 
        costArgs=list(...), envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


#' Resampling-based prediction error for sparse LTS regression models
#' 
#' Estimate the prediction error of a previously fit sparse least trimmed 
#' squares regression model via (repeated) 
#' \eqn{K}-fold cross-validation, (repeated) random splitting (also known as 
#' random subsampling or Monte Carlo cross-validation), or the bootstrap.  If 
#' the model fit contains estimates over a grid of values for the penalty 
#' parameter (i.e., for objects of class \code{"sparseLTSGrid"}), then in each 
#' iteration the optimal model is selected from the training data and used to 
#' make predictions for the test data.
#' 
#' @method perry sparseLTS
#' 
#' @param object  the model fit for which to estimate the prediction error.
#' @param splits  an object of class \code{"cvFolds"} (as returned by 
#' \code{\link[perry]{cvFolds}}) or a control object of class 
#' \code{"foldControl"} (see \code{\link[perry]{foldControl}}) defining the 
#' folds of the data for (repeated) \eqn{K}-fold cross-validation, an object of 
#' class \code{"randomSplits"} (as returned by 
#' \code{\link[perry]{randomSplits}}) or a control object of class 
#' \code{"splitControl"} (see \code{\link[perry]{splitControl}}) defining 
#' random data splits, or an object of class \code{"bootSamples"} (as returned 
#' by \code{\link[perry]{bootSamples}}) or a control object of class 
#' \code{"bootControl"} (see \code{\link[perry]{bootControl}}) defining 
#' bootstrap samples.
#' @param fit  a character string specifying for which fit to estimate the 
#' prediction error.  Possible values are \code{"reweighted"} (the default) for 
#' the prediction error of the reweighted fit, \code{"raw"} for the prediction 
#' error of the raw fit, or \code{"both"} for the prediction error of both 
#' fits.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' vectors to be passed as its first two arguments, the first corresponding to 
#' the observed values of the response and the second to the predicted values, 
#' and must return a non-negative scalar value.  The default is to use the root 
#' trimmed mean squared prediction error for robust models (see 
#' \code{\link[perry]{cost}}).
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed to the prediction loss 
#' function \code{cost}.
#' 
#' @returnClass perry 
#' @returnItem pe  a numeric vector containing the estimated prediction errors 
#' for the requested model fits.  In case of more than one replication, this 
#' gives the average value over all replications.
#' @returnItem se  a numeric vector containing the estimated standard errors of 
#' the prediction loss for the requested model fits.
#' @returnItem reps  a numeric matrix in which each column contains the 
#' estimated prediction errors from all replications for the requested model 
#' fits.  This is only returned in case of more than one replication.
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error.
#' @returnItem y  the response.
#' @returnItem yHat  a list containing the predicted values from all 
#' replications.
#' @returnItem seed  the seed of the random number generator before estimation 
#' of the prediction error.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{sparseLTS}}, \code{\link{sparseLTSGrid}}, 
#' \code{\link{predict.sparseLTS}}, \code{\link[perry]{perry}}, 
#' \code{\link[perry]{cost}}
#' 
#' @example inst/doc/examples/example-perry.sparseLTS.R
#' 
#' @keywords utilities robust
#' 
#' @export
#' @import perry

perry.sparseLTS <- function(object, splits = foldControl(), 
        fit = c("reweighted", "raw", "both"), cost = rtmspe, 
        seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perry")
    fit <- match.arg(fit)
    if(is.null(x <- object$x) || is.null(y <- object$y)) {
        if(is.null(x)) x <- try(model.matrix(object$terms), silent=TRUE)
        if(is.null(y)) y <- try(model.response(object$terms), silent=TRUE)
        if(inherits(x, "try-error") || inherits(y, "try-error")) {
            stop("model data not available")
        }
    }
    # predictor matrix is stored with column for intercept (if any)
    x <- removeIntercept(x)
    ## prepare prediction error estimation
    # extract function call for fitting the model
    call <- object$call
    # if the model was fitted with formula method, 'formula' and 'data' 
    # arguments are removed from call and 'x' and 'y' are used instead
    call$formula <- NULL
    call$data <- NULL
    call$intercept <- object$intercept
    ## call function perryFit() to perform prediction error estimation
    out <- perryFit(call, x=x, y=y, splits=splits, predictArgs=list(fit=fit), 
        cost=cost, costArgs=list(...), envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


#' @rdname perry.sparseLTS
#' @method perry sparseLTSGrid
#' @export

perry.sparseLTSGrid <- function(object, splits = foldControl(), 
        fit = c("reweighted", "raw", "both"), cost = rtmspe, 
        seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perry")
    fit <- match.arg(fit)
    ## call "sparseLTS" method and adjust prediction error names if necessary
    out <- perry.sparseLTS(object=object, splits=splits, 
        fit=fit, cost=cost, seed=seed, ...)
    if(fit == "both") peNames(out) <- c("reweighted", "raw")
    out$call <- matchedCall
    out
}
