# ---------------------------------------
# Authors: Andreas Alfons
#          Erasmus Universiteit Rotterdam
#
#          Viktoria Oellerer
#          KU Leuven
# ---------------------------------------

#' Sparse MM-estimator of regression
#'
#' Compute the MM-estimator of regression with an \eqn{L_{1}}{L1} penalty on
#' the regression coefficients, which allows for sparse model estimates.
#'
#' @aliases print.sparseMM
#'
#' @encoding utf8
#'
#' @param formula  a formula describing the model.
#' @param data  an optional data frame, list or environment (or object coercible
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in
#' the model.  If not found in data, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{sparseMM} is called.
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a numeric vector of non-negative values to be used as penalty
#' parameter.  Defaults to \code{NULL}, in which case the optimal value of the
#' penalty parameter according to BIC is determined in each weighted lasso fit
#' in the iterative algorithm.
#' @param lambdaS the values of the penalty parameter for in initial sparse
#' S-estimator can be specified separately. This should be a numeric vector of
#' the same length as \code{lambda}.
#' @param standardize  a logical indicating whether the response and predictor
#' variables should be robustly standardized (the default is \code{TRUE}). See
#' \code{\link[=standardize]{robStandardize}}.
#' @param tuning.psi  numeric; tuning constant for Tukey's bisquare function
#' that determines the efficiency. The default value yields 95\% efficiency for
#' the unpenalized MM-estimator at the regression model with normal errors.
#' Use value 3.443689 for 85\% efficiency.
#' @param max.it  integer; the maximum number of iteratively reweighted lasso
#' iterations.
#' @param nsamp,nistep,tuning.chi,bb,nfpi  see \code{\link{sparseS}}.
#' @param tol  a small positive numeric value giving the tolerance for
#' convergence.
#' @param eps  a small positive numeric value used to determine whether the
#' variability within a variable is too small (an effective zero).
#' @param use.Gram  a logical indicating whether the Gram matrix of the
#' explanatory variables should be precomputed in the weighted lasso fits. If
#' the number of variables is large, computation may be faster when this is set
#' to \code{FALSE}.  The default is to use \code{TRUE} if the number of
#' variables is smaller than the number of observations and smaller than 100,
#' and \code{FALSE} otherwise.
#' @param scaleType  a character string specifying the type of residual scale
#' estimate.  Possible values are \code{"S"} (the default; the scale estimate
#' from the initial S-estimator), \code{"M"} (re-estimated M-scale of the final
#' residuals) and \code{"tau"} (\eqn{\tau}{tau}-estimate of scale of the final
#' residuals, see \code{\link[robustbase]{scaleTau2}}).
#' @param crit  a character string specifying the optimality criterion to be
#' used for selecting the final model.  Currently only \code{"BIC"} for
#' the Bayes information criterion is implemented.
#' @param ncores  a positive integer giving the number of processor cores to be
#' used for parallel computing (the default is 1 for no parallelization).  If
#' this is set to \code{NA}, all available processor cores are used.  Parallel
#' computing is implemented on the C++ level via OpenMP
#' (\url{http://openmp.org/}).
#' @param seed  optional initial seed for the random number generator (see
#' \code{\link{.Random.seed}}).
#' @param model  a logical indicating whether the data \code{x} and \code{y}
#' should be added to the return object.  A column of ones is added to \code{x}
#' to account for the intercept.
#' @param \dots  additional arguments to be passed down.
#'
#' @return
#' An object of class \code{"sparseMM"} (inheriting from class
#' \code{"penModel"}) with the following components:
#' \item{lambda}{a numeric vector giving the values of the penalty
#' parameter.}
#' \item{rweights}{an integer vector or matrix containing the respective
#' robustness weights.}
#' \item{objective}{a numeric vector giving the respective values of the
#' sparse MM objective function, i.e., the \eqn{L_{1}}{L1} penalized sums of
#' the M-loss of the scaled residuals.}
#' \item{coefficients}{a numeric vector or matrix containing the
#' respective coefficient estimates.}
#' \item{fitted.values}{a numeric vector or matrix containing the
#' respective fitted values of the response.}
#' \item{residuals}{a numeric vector or matrix containing the
#' respective residuals.}
#' \item{scale}{a numeric vector giving the robust scale estimates of the
#' corresponding residuals.}
#' \item{df}{an integer vector giving the respective degrees of freedom
#' of the obtained model fits, i.e., the number of nonzero coefficient
#' estimates.}
#' \item{standardize}{a logical indicating whether the predictor
#' variables were robustly standardized.}
#' \item{init.S}{a list containing the results of the initial sparse
#' S-estimator.}
#' \item{crit}{an object of class \code{"bicSelect"} containing the BIC
#' values and indicating the final model (only returned if argument \code{crit}
#' is \code{"BIC"} and argument \code{lambda} contains more than one value for
#' the penalty parameter).}
#' \item{muX}{a numeric vector containing the center estimates of the
#' predictors if \code{standardize=TRUE} and a numeric vectors of zeros
#' otherwise.}
#' \item{sigmaX}{a numeric vector containing the scale estimates of the
#' predictors if \code{standardize=TRUE} and a numeric vectors of ones
#' otherwise.}
#' \item{muY}{numeric; the center estimate of the response if
#' \code{standardize=TRUE} and a numeric vectors of zeros otherwise.}
#' \item{sigmaY}{numeric; the scale estimate of the response if
#' \code{standardize=TRUE} and a numeric vectors of ones otherwise.}
#' \item{x}{the predictor matrix (if \code{model} is \code{TRUE}).}
#' \item{y}{the response variable (if \code{model} is \code{TRUE}).}
#' \item{call}{the matched function call.}
#'
#' @author Andreas Alfons and Viktoria \enc{Öllerer}{Oellerer}
#'
#' @seealso \code{\link[=coef.penModel]{coef}},
#' \code{\link[=fitted.penModel]{fitted}},
#' \code{\link[=plot.penModel]{plot}},
#' \code{\link[=predict.penModel]{predict}},
#' \code{\link[=residuals.penModel]{residuals}},
#' \code{\link[=weights.penModel]{weights}},
#' \code{\link[robustbase]{lmrob}}
#'
#' @keywords regression robust
#'
#' @export
#' @import parallel

sparseMM <- function(x, ...) UseMethod("sparseMM")


#' @rdname sparseMM
#' @method sparseMM formula
#' @export

sparseMM.formula <- function(formula, data, ...) {
  ## get function call
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("sparseMM")
  ## prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if(is.empty.model(mt)) stop("empty model")
  ## extract response and candidate predictors from model frame
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  ## call default method
  fit <- sparseMM(x, y, ...)
  # return results
  fit$call <- matchedCall  # add call to return object
  fit$terms <- mt          # add model terms to return object
  fit
}


#' @rdname sparseMM
#' @method sparseMM default
#' @export
sparseMM.default <- function(x, y, lambdaMin = c(0, 0), sMax = NA,
                             standardize = TRUE, tuning.psi = 4.685061,
                             max.it = 500, nsamp = c(500, 10), nistep = 10,
                             tuning.chi = 1.547645, bb = NULL,
                             nfpi = c(10, 500), tol = .Machine$double.eps^0.5,
                             eps = .Machine$double.eps, use.Gram,
                             scaleType = c("S", "M", "tau"), crit = "BIC",
                             ncores = 1, seed = NULL, model = TRUE, ...) {

  ## initializations
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("sparseMM")
  n <- length(y)
  x <- addColnames(as.matrix(x))
  d <- dim(x)
  if(!isTRUE(n == d[1])) stop(sprintf("'x' must have %d rows", n))
  # TODO: sensible default for lambdaMin as small fraction of lambda0() if n < p
  # TODO: use function lambdaToLambdaS()
  lambdaMin <- rep(lambdaMin, length.out=2)
  if(!is.numeric(lambdaMin) || any(!is.finite(lambdaMin))) {
    lambdaMin <- formals()$lambdaMin
    warning("missing or infinite value of 'lambdaMin'; using default value")
  }
  if(is.na(sMax)) sMax <- d[2]  # suitable value for n < p determined in C++
  if(!is.numeric(sMax) || is.infinite(sMax) || sMax < 0 || sMax > d[2]) {
    sMax <- d[2]  # use default value
    warning("invalid value of 'sMax'; using default value")
  } else sMax <- as.integer(sMax)
  standardize <- isTRUE(standardize)
  tuning.psi <- rep(tuning.psi, length.out=1)
  if(!is.numeric(tuning.psi) || !is.finite(tuning.psi)) {
    tuning.psi <- formals()$tuning.psi
    warning("missing or infinite value of 'tuning.psi'; using default value")
  }
  max.it <- rep(max.it, length.out=1)
  if(!is.numeric(max.it) || !is.finite(max.it)) {
    max.it <- formals()$max.it
    warning("missing or infinite value of 'max.it'; using default value")
  } else max.it <- as.integer(max.it)
  nsamp <- rep(nsamp, length.out=2)
  if(!is.numeric(nsamp) || any(!is.finite(nsamp))) {
    nsamp <- formals()$nsamp
    warning("missing or infinite values in 'nsamp'; using default values")
  } else nsamp <- as.integer(nsamp)
  nistep <- rep(as.integer(nistep), length.out=1)
  if(!is.numeric(nistep) || !is.finite(nistep)) {
    nistep <- formals()$nistep
    warning("missing or infinite value of 'nistep'; using default value")
  } else nistep <- as.integer(nistep)
  tuning.chi <- rep(tuning.chi, length.out=1)
  if(!is.numeric(tuning.chi) || !is.finite(tuning.chi)) {
    tuning.chi <- formals()$tuning.chi
    warning("missing or infinite value of 'tuning.chi'; using default value")
  }
  bbDefault <- Tbsb(tuning.chi, 1)/tuning.chi * Mrho(tuning.chi, tuning.chi)
  if(is.null(bb)) bb <- bbDefault
  bb <- rep(bb, length.out=1)
  if(!is.numeric(bb) || !is.finite(bb)) {
    bb <- bbDefault
    warning("missing or infinite value of 'bb'; using default value")
  }
  nfpi <- rep(nfpi, length.out=2)
  if(!is.numeric(nfpi) || any(!is.finite(nfpi))) {
    nfpi <- formals()$nfpi
    warning("missing or infinite values in 'nfpi'; using default values")
  } else nfpi <- as.integer(nfpi)
  tol <- rep(tol, length.out=1)
  if(!is.numeric(tol) || !is.finite(tol)) {
    tol <- formals()$tol
    warning("missing or infinite value of 'tol'; using default value")
  }
  eps <- rep(eps, length.out=1)
  if(!is.numeric(eps) || !is.finite(eps)) {
    eps <- formals()$eps
    warning("missing or infinite value of 'eps'; using default value")
  }
  if(missing(use.Gram)) use.Gram <- d[2] >= min(n, 100)
  else use.Gram <- isTRUE(use.Gram)
  scaleType <- match.arg(scaleType)
  crit <- "none"  # CV not yet implemented
  ncores <- rep(ncores, length.out=1)
  if(is.na(ncores)) ncores <- detectCores()  # use all available cores
  if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
    ncores <- 1  # use default value
    warning("invalid value of 'ncores'; using default value")
  } else ncores <- as.integer(ncores)

  ## obtain random subsets with only 3 observations
  if(!is.null(seed)) set.seed(seed)  # set seed of random number generator
  subsets <- randomSubsets(n, 3, nsamp[1])

  # robustly standardize data
  if(standardize) {
    z <- robStandardize(y, eps=eps, ...)   # standardize response
    xs <- robStandardize(x, eps=eps, ...)  # standardize predictors
    muY <- attr(z, "center")
    sigmaY <- attr(z, "scale")
    muX <- attr(xs, "center")
    sigmaX <- attr(xs, "scale")
  } else {
    z <- y
    xs <- x
    muY <- 0
    sigmaY <- 1
    muX <- rep.int(0, d[2])
    sigmaX <- rep.int(1, d[2])
  }

  ## call internal function
  fit <- fastSparseMM(x=xs, y=z, lambdaMin=lambdaMin, sMax=sMax,
                      tuning.psi=tuning.psi, max.it=max.it, initial=subsets,
                      nistep=nistep, nkeep=nsamp[2], tuning.chi=tuning.chi,
                      bb=bb, nfpi=nfpi, tol=tol, eps=eps, use.Gram=use.Gram,
                      ncores=ncores)

  ## construct return object
  x <- addIntercept(x)
  # backtransform coefficients, residuals and scale
  coefficients <- backtransform(fit$coefficients, centerY=muY, scaleY=sigmaY,
                                centerX=muX, scaleX=sigmaX)
  residuals <- fit$residuals * sigmaY
  if(scaleType == "M") {
    scale <- Mscale(residuals, k=tuning.chi, b=bb, nfpi=nfpi[2], tol=tol)
  } else if(scaleType == "tau") scale <- scaleTau2(residuals)
  else scale <- fit$scaleS*sigmaY
  df <- modelDf(coefficients, tol)
  out <- list(lambda=fit$lambda,
              rweights=copyNames(from=y, to=fit$weights),
              objective=fit$objective,
              coefficients=copyNames(from=x, to=coefficients),
              fitted.values=copyNames(from=y, to=y-residuals),
              residuals=copyNames(from=y, to=residuals),
              scale=scale, df=df, standardize=standardize)
  # add information on initial S-estimator
  coefficientsS <- backtransform(fit$coefficientsS,
                                 centerY=muY, scaleY=sigmaY,
                                 centerX=muX, scaleX=sigmaX)
  residualsS <- fit$residualsS * sigmaY
  dfS <- modelDf(coefficientsS, tol)
  out$init.S <- list(lambda=fit$lambdaS,
                     rweights=copyNames(from=y, to=fit$weightsS),
                     objective=fit$objectiveS,
                     coefficients=copyNames(from=x, to=coefficientsS),
                     fitted.values=copyNames(from=y, to=y-residualsS),
                     residuals=copyNames(from=y, to=residualsS),
                     scale=fit$scaleS*sigmaY, df=dfS)

  ## assign class
  class(out) <- c("sparseMM", "penModel")

  ## add information on the optimal model
  if(crit == "BIC") out$crit <- bicSelect(out)

  ## add center and scale estimates
  out[c("muX", "sigmaX", "muY", "sigmaY")] <- list(muX, sigmaX, muY, sigmaY)
  ## add model data if requested
  if(isTRUE(model)) out[c("x", "y")] <- list(x=x, y=y)

  ## return results
  out$call <- matchedCall
  out
}


## internal function to compute sparse MM-estimator
fastSparseMM <- function(x, y, lambdaMin, sMax, tuning.psi = 4.685061,
                         max.it = 500, initial, nistep = 10, nkeep = 10,
                         tuning.chi = 1.547645, bb = 0.5, nfpi = c(10, 500),
                         tol = .Machine$double.eps^0.5,
                         eps = .Machine$double.eps, use.Gram = TRUE,
                         ncores = 1, drop = TRUE) {
  # call C++ function
  fit <- .Call("R_fastSparseMM", R_x=x, R_y=y, R_lambdaMinM=lambdaMin[1],
               R_lambdaMinS=lambdaMin[2], R_sMax=sMax,
               R_k=tuning.psi, R_rMax=max.it, R_initial=initial,
               R_nistep=nistep, R_nkeep=nkeep, R_kS=tuning.chi, R_b=bb,
               R_nfpi=nfpi[1], R_nfpiMax=nfpi[2], R_tol=tol, R_eps=eps,
               R_useGram=use.Gram, R_ncores=ncores, PACKAGE="robustHD")
  if(drop) {
    # drop the dimension of selected components
    which <- c("weights", "residuals", "weightsS", "residualsS")
    fit[which] <- lapply(fit[which], drop)
  }
  # return MM-estimate
  fit
}


## internal function to compute value of lambdaS such that lambda and lambdaS
## achieve the same amount of shrinkage
lambdaToLambdaS <- function(lambda, k = 4.685061, kS = 1.547645) {
  # 2*lambda*E[psi1(Z)/Z] / (E[psi1(Z)*Z]*E[psi2(Z)/Z])
  2*lambda/(2/kS^4*ksiint(kS, 3, 1)-4/kS^2*ksiint(kS, 2, 1)+2*ksiint(kS, 1, 1))*
    (2/kS^4*ksiint(kS, 2, 1)-4/kS^2*ksiint(kS, 1, 1)+2*(2*pnorm(kS)-1)) /
    (2/k^4*ksiint(k, 2, 1)-4/k^2*ksiint(k, 1, 1)+2*(2*pnorm(k)-1))
}
