# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Predict from a sequence of regression models
#' 
#' Make predictions from a sequence of regression models such as submodels 
#' along a robust or groupwise least angle regression sequence.  For 
#' autoregressive time series models with exogenous inputs, \eqn{h}-step ahead 
#' forecasts are performed.
#' 
#' For \code{predict.seqModel}, the \code{newdata} argument defaults to the 
#' matrix of predictors used to fit the model such that the fitted values are 
#' computed.
#' 
#' For autoregressive time series models with exogenous inputs with forecast 
#' horizon \eqn{h}, the \eqn{h} most recent observations of the predictors are 
#' omitted from fitting the model since there are no corresponding values for 
#' the response.  Hence the \code{newdata} argument for \code{predict.tslarsP} 
#' and \code{predict.tslars} defaults to those \eqn{h} observations of the 
#' predictors.
#' 
#' @method predict seqModel
#' @aliases predict.rlars predict.grplars
#' 
#' @param object  the model fit from which to make predictions.
#' @param newdata  new data for the predictors.  If the model fit was computed 
#' with the formula method, this should be a data frame from which to extract 
#' the predictor variables.  Otherwise this should be a matrix containing the 
#' same variables as the predictor matrix used to fit the model (including a 
#' column of ones to account for the intercept).
#' @param p  an integer giving the lag length for which to make predictions 
#' (the default is to use the optimal lag length). 
#' @param s  an integer vector giving the steps of the submodels for which to 
#' make predictions (the default is to use the optimal submodel).
#' @param \dots  for the \code{"tslars"} method, additional arguments to be 
#' passed down to the \code{"tslarsP"} method.  For the other methods, 
#' additional arguments to be passed down to the respective method of 
#' \code{\link[=coef.seqModel]{coef}}.
#' 
#' @return  
#' If only one submodel is requested, a numeric vector containing the 
#' corresponding predicted values.
#' 
#' If multiple submodels are requested, a numeric matrix in which each column 
#' contains the predicted values from the corresponding submodel.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{predict}}, \code{\link{rlars}}, 
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}}, 
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}}
#' 
#' @example inst/doc/examples/example-predict.rlars.R
#' 
#' @keywords regression
#' 
#' @export

predict.seqModel <- function(object, newdata, s = NA, ...) {
    ## initializations
    coef <- coef(object, s=s, ...)  # extract coefficients
    # extract or check new data
    d <- dim(coef)
    terms <- delete.response(object$terms)  # extract terms for model matrix
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(newdata <- object$x)) {
            newdata <- try(model.matrix(terms), silent=TRUE)
            if(inherits(newdata, "try-error")) stop("model data not available")
        }
    } else {
        # interpret vector as row
        if(is.null(dim(newdata))) newdata <- t(newdata)
        # check dimensions if model was not specified with a formula, 
        # otherwise use the terms object to extract model matrix
        if(is.null(terms)) {
            newdata <- as.matrix(newdata)
            # add a column of ones to the new data matrix 
            # (unless it already contains intercept column)
            newdata <- addIntercept(newdata, check=TRUE)
            # check dimensions of new data
            p <- if(is.null(d)) length(coef) else d[1]
            if(ncol(newdata) != p) {
                stop(sprintf("new data must have %d columns", p))
            }
        } else newdata <- model.matrix(terms, as.data.frame(newdata))
    }
    ## compute predictions
    # ensure that a vector is returned if only one fit is requested
    out <- newdata %*% coef
    if(is.null(d)) out <- drop(out)
    out
}


#' @rdname predict.seqModel
#' @method predict optSeqModel
#' @export

predict.optSeqModel <- function(object, newdata, ...) {
    ## initializations
    coef <- coef(object, ...)  # extract coefficients
    # extract or check new data
    terms <- delete.response(object$terms)  # extract terms for model matrix
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(newdata <- object$x)) {
            newdata <- try(model.matrix(terms), silent=TRUE)
            if(inherits(newdata, "try-error")) stop("model data not available")
        }
    } else {
        # interpret vector as row
        if(is.null(dim(newdata))) newdata <- t(newdata)
        # check dimensions if model was not specified with a formula, 
        # otherwise use the terms object to extract model matrix
        if(is.null(terms)) {
            newdata <- as.matrix(newdata)
            # add a column of ones to the new data matrix 
            # (unless it already contains intercept column)
            newdata <- addIntercept(newdata, check=TRUE)
            # check dimensions of new data
            p <- length(coef)
            if(ncol(newdata) != p) {
                stop(sprintf("new data must have %d columns", p))
            }
        } else newdata <- model.matrix(terms, as.data.frame(newdata))
    }
    ## compute predictions (ensure that a vector is returned)
    drop(newdata %*% coef)
}


#' @rdname predict.seqModel
#' @method predict tslarsP
#' @export

predict.tslarsP <- function(object, newdata, ...) {
    ## initializations
    coef <- coef(object, ...)  # extract coefficients
    d <- dim(coef)
    terms <- object$terms  # extract terms for model matrix
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(x <- object$x) || is.null(y <- object$y)) {
            if(is.null(x)) x <- try(model.matrix(object$terms), silent=TRUE)
            if(is.null(y)) y <- try(model.response(object$terms), silent=TRUE)
            if(inherits(x, "try-error") || inherits(y, "try-error")) {
                stop("model data not available")
            }
        }
        newdata <- newdataBlocks(x, y, object$h, object$p)
    } else {
        # interpret vector as row
        if(is.null(dim(newdata))) newdata <- t(newdata)
        # check dimensions if model was not specified with a formula, 
        # otherwise use the terms object to extract model matrix
        if(is.null(terms)) {
            newdata <- as.matrix(newdata)
            # add a column of ones to the new data matrix 
            # (unless it already contains intercept column)
            newdata <- addIntercept(newdata, check=TRUE)
            # check dimensions of new data
            p <- if(is.null(d)) length(coef) else d[1]
            if(ncol(newdata) != p) {
                stop(sprintf("new data must have %d columns", p))
            }
        } else {
            mf <- model.frame(terms, as.data.frame(newdata))
            y <- model.response(mf)
            x <- model.matrix(terms, mf)
            if(attr(terms, "intercept")) x <- x[, -1, drop=FALSE]
            newdata <- tsBlocks(x, y, object$p)
        }
    }
    ## compute predictions
    # ensure that a vector is returned if only one fit is requested
    out <- newdata %*% coef
    if(is.null(d)) out <- drop(out)
    out
    
}


#' @rdname predict.seqModel
#' @method predict tslars
#' @export

predict.tslars <- function(object, newdata, p, ...) {
    ## initializations
    # check lag length
    if(missing(p) || !is.numeric(p) || length(p) == 0) {
        p <- object$pOpt
    } else p <- p[1]
    pMax <- object$pMax
    if(p < 1) {
        p <- 1
        warning("lag length too small, using lag length 1")
    } else if(p > pMax) {
        p <- pMax
        warning(sprintf("lag length too large, using maximum lag length %d", p))
    }
    ## compute predictions
    # if missing, construct newdata
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(x <- object$x) || is.null(y <- object$y)) {
            if(is.null(x)) x <- try(model.matrix(object$terms), silent=TRUE)
            if(is.null(y)) y <- try(model.response(object$terms), silent=TRUE)
            if(inherits(x, "try-error") || inherits(y, "try-error")) {
                stop("model data not available")
            }
        }
        # extract model for specified lag length and add original data
        object <- object$pFit[[p]]
        object$x <- x
        object$y <- y
        predict(object, ...)
    } else predict(object$pFit[[p]], newdata, ...)
}


#' Predict from sparse LTS regression models
#' 
#' Make predictions from sparse least trimmed squares regression models.
#' 
#' The \code{newdata} argument defaults to the matrix of predictors used to fit 
#' the model such that the fitted values are computed.
#' 
#' @method predict sparseLTS
#' 
#' @param object  the model fit from which to make predictions.
#' @param newdata  new data for the predictors.  If the model fit was computed 
#' with the formula method, this should be a data frame from which to extract 
#' the predictor variables.  Otherwise this should be a matrix containing the 
#' same variables as the predictor matrix used to fit the model (possibly 
#' including a column of ones in case of a model with intercept).
#' @param s  an integer vector giving the indices of the models for which to 
#' make predictions.  If \code{fit} is \code{"both"}, this can be a list 
#' with two components, with the first component giving the indices of the 
#' reweighted fits and the second the indices of the raw fits.  The default is 
#' to use the optimal model for each of the requested estimators.  Note that 
#' the optimal models may not correspond to the same value of the penalty 
#' parameter for the reweighted and the raw estimator.
#' @param fit  a character string specifying for which fit to make 
#' predictions.  Possible values are \code{"reweighted"} (the default) for 
#' predicting values from the reweighted fit, \code{"raw"} for predicting 
#' values from the raw fit, or \code{"both"} for predicting values from both 
#' fits.
#' @param \dots  currently ignored.
#' 
#' @return  
#' If predictions for only one model are requested, they are returned in the 
#' form of a numeric vector.
#' 
#' Otherwise a numeric matrix is returned in which each column contains the 
#' predicted values from the corresponding model.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[stats]{predict}}, \code{\link{sparseLTS}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-predict.sparseLTS.R
#' 
#' @keywords regression
#' 
#' @export

predict.sparseLTS <- function(object, newdata, 
        fit = c("reweighted", "raw", "both"), ...) {
    ## initializations
    coef <- coef(object, fit=fit)  # extract coefficients
    d <- dim(coef)
    terms <- delete.response(object$terms)  # extract terms for model matrix
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(newdata <- object$x)) {
            newdata <- try(model.matrix(terms), silent=TRUE)
            if(inherits(newdata, "try-error")) stop("model data not available")
        }
    } else {
        # interpret vector as row
        if(is.null(dim(newdata))) newdata <- t(newdata)
        # check dimensions if model was not specified with a formula, 
        # otherwise use the terms object to extract model matrix
        if(is.null(terms)) {
            newdata <- as.matrix(newdata)
            if(object$intercept) {
                # if model has an intercept, add a column of ones to the new 
                # data matrix (unless it already contains intercept column)
                newdata <- addIntercept(newdata, check=TRUE)
            }
            # check dimensions of new data
            p <- if(is.null(d)) length(coef) else d[1]
            if(ncol(newdata) != p) {
                stop(sprintf("new data must have %d columns", p))
            }
        } else newdata <- model.matrix(terms, as.data.frame(newdata))
    }
    ## compute predictions
    # ensure that a vector is returned if only one fit is requested
    out <- newdata %*% coef
    if(is.null(d)) out <- drop(out)
    out
}


#' @rdname predict.sparseLTS
#' @method predict sparseLTSGrid
#' @export

predict.sparseLTSGrid <- function(object, newdata, s = NA, 
        fit = c("reweighted", "raw", "both"), ...) {
    ## initializations
    coef <- coef(object, s=s, fit=fit)  # extract coefficients
    d <- dim(coef)
    terms <- delete.response(object$terms)  # extract terms for model matrix
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(newdata <- object$x)) {
            newdata <- try(model.matrix(terms), silent=TRUE)
            if(inherits(newdata, "try-error")) stop("model data not available")
        }
    } else {
        # interpret vector as row
        if(is.null(dim(newdata))) newdata <- t(newdata)
        # check dimensions if model was not specified with a formula, 
        # otherwise use the terms object to extract model matrix
        if(is.null(terms)) {
            newdata <- as.matrix(newdata)
            # add a column of ones to the new data matrix 
            # (unless it already contains intercept column)
            newdata <- addIntercept(newdata, check=TRUE)
            # check dimensions of new data
            p <- if(is.null(d)) length(coef) else d[1]
            if(ncol(newdata) != p) {
                stop(sprintf("new data must have %d columns", p))
            }
        } else newdata <- model.matrix(terms, as.data.frame(newdata))
    }
    ## compute predictions
    # ensure that a vector is returned if only one fit is requested
    out <- newdata %*% coef
    if(is.null(d)) out <- drop(out)
    out
}


#' @rdname predict.sparseLTS
#' @method predict optSparseLTSGrid
#' @export

predict.optSparseLTSGrid <- predict.sparseLTS
