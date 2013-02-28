# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

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
#' @seealso \code{\link[stats]{predict}}, \code{\link{sparseLTS}}
#' 
#' @example inst/doc/examples/example-predict.sparseLTS.R
#' 
#' @keywords regression
#' 
#' @export

predict.sparseLTS <- function(object, newdata, s = NA, 
                              fit = c("reweighted", "raw", "both"), 
                              ...) {
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
