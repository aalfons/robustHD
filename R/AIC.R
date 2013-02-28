# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## class 'lmrob'

AIC.lmrob <- AIC.lmrob.S <- function(object, ..., k = 2) {
    n <- length(residuals(object))  # number of observations
    # compute AIC with the same terms as R does for linear models
    n * (log(2 * pi) + 1 + log(object$scale^2)) + (object$rank) * k
}

BIC.lmrob <- BIC.lmrob.S <- function(object, ...) {
    n <- length(residuals(object))  # number of observations
    AIC(object, ..., k=log(n))      # call AIC method with penalty for BIC
}


## class 'rlm'

AIC.rlm <- function(object, ..., k = 2) {
    n <- length(residuals(object))  # number of observations
    # compute AIC with the same terms as R does for linear models
    n * (log(2 * pi) + 1 + log(object$s^2)) + (object$rank) * k
}

BIC.rlm <- function(object, ...) {
    n <- length(residuals(object))  # number of observations
    AIC(object, ..., k=log(n))      # call AIC method with penalty for BIC
}


#' Information criteria for sparse LTS regression models
#' 
#' Compute the Akaike or Bayes information criterion for sparse least trimmed 
#' squares regression models based on the robust residual scale estimate.
#' 
#' The information criteria are computed as
#' \eqn{n (\log(2 \pi) + 1 + \log(\hat{\sigma}^2)) + df k}{n (log(2 pi) + 1 + log(sigma^2)) + df k}, 
#' where \eqn{n} denotes the number of observations, \eqn{\hat{\sigma}}{sigma} 
#' is the robust residual scale estimate, \eqn{df} is the number of nonzero 
#' coefficient estimates, and \eqn{k} is penalty per parameter.  The usual 
#' definition of the AIC uses \eqn{k = 2}, whereas the BIC uses 
#' \eqn{k = \log(n)}{k = log(n)}.  Consequently, the former is used as the 
#' default penalty of the \code{AIC} method, whereas the \code{BIC} method calls 
#' the \code{AIC} method with the latter penalty.
#' 
#' @method AIC sparseLTS
#' 
#' @param object  the model fit for which to compute the information criterion.
#' @param \dots  for the \code{BIC} method, additional arguments to be passed 
#' down to the \code{AIC} method.  For the \code{AIC} method, additional 
#' arguments are currently ignored.
#' @param fit  a character string specifying for which fit to compute the 
#' information criterion.  Possible values are \code{"reweighted"} (the 
#' default) for the information criterion of the reweighted fit, \code{"raw"} 
#' for the information criterion of the raw fit, or \code{"both"} for the 
#' information criteria of both fits.
#' @param k  a numeric value giving the penalty per parameter to be used.  The 
#' default is to use \eqn{2} as in the classical definition of the AIC.
#' 
#' @return  
#' A numeric vector or matrix giving the information criteria for the requested 
#' fits.
#' 
#' @note Computing information criteria for several objects supplied via the 
#' \code{\dots} argument (as for the default methods of \code{\link[stats]{AIC}} 
#' and \code{BIC}) is currently not implemented.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Akaike, H. (1970) Statistical predictor identification. \emph{Annals of the 
#' Institute of Statistical Mathematics}, \bold{22}(2), 203--217.
#' 
#' Schwarz, G. (1978) Estimating the dimension of a model. \emph{The Annals of 
#' Statistics}, \bold{6}(2), 461--464.
#' 
#' @seealso \code{\link[stats]{AIC}}, \code{\link{sparseLTS}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-AIC.R
#' 
#' @keywords regression
#' 
#' @import stats
#' @export

AIC.sparseLTS <- function(object, ..., 
                          fit = c("reweighted", "raw", "both"), 
                          k = 2) {
  n <- length(residuals(object))  # number of observations
  fit <- match.arg(fit)
  if(fit == "reweighted") {
    sigma2 <- object$scale^2
    df <- object$df
  } else if(fit == "raw") {
    sigma2 <- object$raw.scale^2
    df <- object$raw.df
  } else {
    sigma2 <- cbind(reweighted=object$scale, raw=object$raw.scale)^2
    df <- cbind(reweighted=object$df, raw=object$raw.df)
  }
  # compute AIC with the same terms as R does for linear models
  n * (log(2 * pi) + 1 + log(sigma2)) + df * k
}


#' @rdname AIC.sparseLTS
#' @method BIC sparseLTS
#' @export

BIC.sparseLTS <- function(object, ...) {
  n <- length(residuals(object))  # number of observations
  AIC(object, ..., k=log(n))      # call AIC method with penalty for BIC
}


## internal function returning an object of class "BIC" so that the optimal 
## model can be retrieved for a sequence of models
bicSelect <- function(object, ...) {
  values <- BIC(object, ...)
  if(is.null(dim(values))) best <- which.min(unname(values))
  else best <- apply(values, 2, which.min)
  bic <- list(values=values, best=best)
  class(bic) <- "bicSelect"
  bic
}
