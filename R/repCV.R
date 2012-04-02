# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Cross-validation for a sequential regression model
#' 
#' Perform (repeated) \eqn{K}-fold cross-validation to estimate the prediction 
#' error of a previously fit sequential regression model such as a robust least 
#' angle regression model.  In each iteration of cross-validation, the optimal 
#' model is thereby selected from the training data and used to make 
#' predictions for the test data.
#' 
#' @method repCV seqModel
#' @aliases repCV.rlars
#' 
#' @param object  the model fit for which to estimate the prediction error.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' vectors to be passed as its first two arguments, the first corresponding to 
#' the observed values of the response and the second to the predicted values, 
#' and must return a non-negative scalar value.  The default is to use the root 
#' mean squared prediction error for non-robust models and the root trimmed 
#' mean squared prediction error for robust models (see 
#' \code{\link[cvTools]{cost}}).
#' @param K  an integer giving the number of groups into which the data should 
#' be split (the default is five).  Keep in mind that this should be chosen 
#' such that all groups are of approximately equal size.  Setting \code{K} 
#' equal to \code{n} yields leave-one-out cross-validation.
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.  This is ignored for for leave-one-out 
#' cross-validation and other non-random splits of the data.
#' @param foldType  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' @param folds  an object of class \code{"cvFolds"} giving the folds of the 
#' data for cross-validation (as returned by 
#' \code{\link[cvTools]{cvFolds}}).  If supplied, this is preferred over 
#' \code{K} and \code{R}.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed to the prediction loss 
#' function \code{cost}.
#' 
#' @returnClass cv 
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem cv  a numeric value giving the estimated prediction error.  For 
#' repeated cross-validation, this gives the average value over all 
#' replications.
#' @returnItem se  a numeric value giving the estimated standard error of 
#' the prediction loss.
#' @returnItem reps  a numeric matrix with one column that contains the 
#' estimated prediction errors from all replications.  This is only returned 
#' for repeated cross-validation.
#' @returnItem seed  the seed of the random number generator before 
#' cross-validation was performed.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{rlars}}, \code{\link{predict.seqModel}}, 
#' \code{\link[cvTools]{cvFolds}}, \code{\link[cvTools]{cost}}
#' 
#' @example inst/doc/examples/example-repCV.rlars.R
#' 
#' @keywords utilities
#' 
#' @export
#' @import cvTools

repCV.seqModel <- function(object, cost, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("repCV")
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
    ## construct call to fit models in cross-validation
    call <- object$call
    # if the model was fitted with formula method, 'formula' and 'data' 
    # arguments are removed from call and 'x' and 'y' are used instead
    call$formula <- NULL
    call$data <- NULL
    ## call function cvFit() to perform cross-validation
    out <- cvFit(call, x=x, y=y, cost=cost, K=K, R=R, foldType=foldType, 
        folds=folds, costArgs=list(...), envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


## workhorse function for classes "sparseLTS" and "sparseLTSGrid"
repCVSparseLTS <- function(object, cost = rtmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), folds = NULL, 
        fit = c("reweighted", "raw", "both"), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("repCV")
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
    ## prepare cross-validation
    # extract function call for fitting the model
    call <- object$call
    # if the model was fitted with formula method, 'formula' and 'data' 
    # arguments are removed from call and 'x' and 'y' are used instead
    call$formula <- NULL
    call$data <- NULL
    call$intercept <- object$intercept
    ## call function cvFit() to perform cross-validation
    out <- cvFit(call, x=x, y=y, cost=cost, K=K, R=R, foldType=foldType, 
        folds=folds, predictArgs=list(fit=fit), costArgs=list(...), 
        envir=parent.frame(), seed=seed)
    if(inherits(object, "sparseLTSGrid") && fit == "both") {
        cvNames(out) <- c("reweighted", "raw")
    }
    out$call <- matchedCall
    out
}


#' Cross-validation for sparse LTS regression models
#' 
#' Estimate the prediction error of a previously fit sparse least trimmed 
#' squares regression model via (repeated) \eqn{K}-fold cross-validation.  If 
#' the model fit contains estimates over a grid of values for the penalty 
#' parameter (i.e., for objects of class \code{"sparseLTSGrid"}), then in each 
#' iteration of cross-validation the optimal model is selected from the 
#' training data and used to make predictions for the test data.
#' 
#' @method repCV sparseLTS
#' 
#' @param object  the model fit for which to estimate the prediction error.
#' @param cost  a robust cost function measuring prediction loss.  It should 
#' expect vectors to be passed as its first two arguments, the first 
#' corresponding to the observed values of the response and the second to the 
#' predicted values, and must return a non-negative scalar value.  The default 
#' is to use the root trimmed mean squared prediction error (see 
#' \code{\link[cvTools]{cost}}).
#' @param K  an integer giving the number of groups into which the data should 
#' be split (the default is five).  Keep in mind that this should be chosen 
#' such that all groups are of approximately equal size.  Setting \code{K} 
#' equal to \code{n} yields leave-one-out cross-validation.
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.  This is ignored for for leave-one-out 
#' cross-validation and other non-random splits of the data.
#' @param foldType  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' @param folds  an object of class \code{"cvFolds"} giving the folds of the 
#' data for cross-validation (as returned by 
#' \code{\link[cvTools]{cvFolds}}).  If supplied, this is preferred over 
#' \code{K} and \code{R}.
#' @param fit  a character string specifying for which fit to estimate the 
#' prediction error.  Possible values are \code{"reweighted"} (the default) for 
#' the prediction error of the reweighted fit, \code{"raw"} for the prediction 
#' error of the raw fit, or \code{"both"} for the prediction error of both 
#' fits.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed to the prediction loss 
#' function \code{cost}.
#' 
#' @returnClass cv 
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem cv  a numeric vector containing the estimated prediction errors 
#' for the requested model fits.  For repeated cross-validation, this contains 
#' the average values over all replications.
#' @returnItem se  a numeric vector containing the estimated standard errors of 
#' the prediction loss for the requested model fits.
#' @returnItem reps  a numeric matrix in which each column contains the 
#' estimated prediction errors from all replications for the requested model 
#' fits.  This is only returned for repeated cross-validation.
#' @returnItem seed  the seed of the random number generator before 
#' cross-validation was performed.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{sparseLTS}}, \code{\link{sparseLTSGrid}}, 
#' \code{\link{predict.sparseLTS}}, \code{\link[cvTools]{cvFolds}}, 
#' \code{\link[cvTools]{cost}}
#' 
#' @example inst/doc/examples/example-repCV.sparseLTS.R
#' 
#' @keywords utilities robust
#' 
#' @export
#' @import cvTools

repCV.sparseLTS <- repCVSparseLTS


#' @rdname repCV.sparseLTS
#' @method repCV sparseLTSGrid
#' @export

repCV.sparseLTSGrid <- repCVSparseLTS
