# ---------------------------------------
# Authors: Andreas Alfons
#          Erasmus Universiteit Rotterdam
#
#          Viktoria Oellerer
#          KU Leuven
# ---------------------------------------

#' Sparse S-estimator of regression
#'
#' Compute the S-estimator of regression with an \eqn{L_{1}}{L1} penalty on
#' the regression coefficients, which allows for sparse model estimates.
#'
#' @aliases print.sparseS
#'
#' @encoding utf8
#'
#' @param formula  a formula describing the model.
#' @param data  an optional data frame, list or environment (or object coercible
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in
#' the model.  If not found in data, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{sparseS} is called.
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a numeric vector of non-negative values to be used as penalty
#' parameter.  Defaults to \code{NULL}, in which case the optimal value of the
#' penalty parameter according to BIC is determined in each weighted lasso fit
#' in the iterative algorithm.
#' @param standardize  a logical indicating whether the response and predictor
#' variables should be robustly standardized (the default is \code{TRUE}).  See
#' \code{\link[=standardize]{robStandardize}}.
#' @param nsamp  a numeric vector giving the number of resamples to be used in
#' the two phases of the algorithm.  The first element gives the number of
#' initial resamples to be used.  The second element gives the number of
#' resamples to keep after the first phase of \code{nistep} I-steps.  For
#' those remaining resamples, additional I-steps are performed until
#' convergence.  The default is to first perform \code{nistep} I-steps on 500
#' initial resamples, and then to keep the 10 resamples with the lowest value
#' of the objective function for additional I-steps until convergence.
#' @param nistep  a positive integer giving the number of I-steps to perform on
#' all resamples in the first phase of the algorithm (the default is to
#' perform 10 I-steps).
#' @param tuning.chi  numeric; tuning constant of Tukey's bisquare function
#' that determines the robustness of the sparse S-estimator.  The default value
#' yields an asymptotic breakdown point of 50\%.  Use value 2.937015 for an
#' asymptotic breakdown point of 25\%.
#' @param bb  numeric; expected value of the M-scale.  The default value is
#' computed under the normal model with tuning constant equal to
#' \code{tuning.chi}.
#' @param nfpi  a numeric vector of length two giving the maximum number of
#' fixed point iterations to be used for computing the M-scale.  The first
#' value is used in the first phase of the algorithm when \code{nistep} I-steps
#' are performed on all resamples, the second value is used in the second phase
#' when I-steps are performed until convergence on the best resamples.
#' @param tol  a small positive numeric value giving the tolerance for
#' convergence.
#' @param eps  a small positive numeric value used to determine whether the
#' variability within a variable is too small (an effective zero).
#' @param use.Gram  a logical indicating whether the Gram matrix of the
#' explanatory variables should be precomputed in the weighted lasso fits.  If
#' the number of variables is large, computation may be faster when this is set
#' to \code{FALSE}.  The default is to use \code{TRUE} if the number of
#' variables is smaller than the number of observations and smaller than 100,
#' and \code{FALSE} otherwise.
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
#' An object of class \code{"sparseS"} (inheriting from class
#' \code{"penModel"}) with the following components:
#' @returnItem lambda  a numeric vector giving the values of the penalty
#' parameter.
#' @returnItem rweights  an integer vector or matrix containing the respective
#' robustness weights.
#' @returnItem objective  a numeric vector giving the respective values of the
#' sparse S objective function, i.e., the \eqn{L_{1}}{L1} penalized sums of
#' the M-scales of the residuals.
#' @returnItem coefficients  a numeric vector or matrix containing the
#' respective coefficient estimates.
#' @returnItem fitted.values  a numeric vector or matrix containing the
#' respective fitted values of the response.
#' @returnItem residuals  a numeric vector or matrix containing the
#' respective residuals.
#' @returnItem scale  a numeric vector giving the minimized M-scale estimates
#' of the corresponding residuals.
#' @returnItem df  an integer vector giving the respective degrees of freedom
#' of the obtained model fits, i.e., the number of nonzero coefficient
#' estimates.
#' @returnItem standardize  a logical indicating whether the predictor
#' variables were robustly standardized.
#' @returnItem crit  an object of class \code{"bicSelect"} containing the BIC
#' values and indicating the final model (only returned if argument \code{crit}
#' is \code{"BIC"} and argument \code{lambda} contains more than one value for
#' the penalty parameter).
#' @returnItem muX  a numeric vector containing the center estimates of the
#' predictors if \code{standardize=TRUE} and a numeric vectors of zeros otherwise.
#' @returnItem sigmaX  a numeric vector containing the scale estimates of the
#' predictors if \code{standardize=TRUE} and a numeric vectors of ones otherwise.
#' @returnItem muY  numeric; the center estimate of the response if
#' \code{standardize=TRUE} and a numeric vectors of zeros otherwise.
#' @returnItem sigmaY  numeric; the scale estimate of the response if
#' \code{standardize=TRUE} and a numeric vectors of ones otherwise.
#' @returnItem x  the predictor matrix (if \code{model} is \code{TRUE}).
#' @returnItem y  the response variable (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#'
#' @author Andreas Alfons and Viktoria \enc{Öllerer}{Oellerer}
#'
#' @seealso \code{\link[=coef.penModel]{coef}},
#' \code{\link[=fitted.penModel]{fitted}},
#' \code{\link[=plot.penModel]{plot}},
#' \code{\link[=predict.penModel]{predict}},
#' \code{\link[=residuals.penModel]{residuals}},
#' \code{\link[=weights.penModel]{weights}}
#'
#' @keywords regression robust
#'
#' @export
#' @import parallel

sparseS <- function(x, ...) UseMethod("sparseS")


#' @rdname sparseS
#' @method sparseS formula
#' @export

sparseS.formula <- function(formula, data, ...) {
  ## get function call
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("sparseS")
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
  fit <- sparseS(x, y, ...)
  # return results
  fit$call <- matchedCall  # add call to return object
  fit$terms <- mt          # add model terms to return object
  fit
}


#' @rdname sparseS
#' @method sparseS default
#' @export
sparseS.default <- function(x, y, lambda = NULL, standardize = TRUE,
                            nsamp = c(500, 10), nistep = 10,
                            tuning.chi = 1.547645, bb = NULL,
                            nfpi = c(10, 500), tol = .Machine$double.eps^0.5,
                            eps = .Machine$double.eps, use.Gram,
                            crit = "BIC", ncores = 1, seed = NULL,
                            model = TRUE, ...) {

  ## initializations
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("sparseS")
  n <- length(y)
  x <- addColnames(as.matrix(x))
  d <- dim(x)
  if(!isTRUE(n == d[1])) stop(sprintf("'x' must have %d rows", n))
  findLambda <- is.null(lambda)
  if((!is.numeric(lambda) || length(lambda) == 0 || any(!is.finite(lambda))) &&
       !findLambda) {
    warning("missing or invalid value of 'lambda'; ",
            "using optimal lambda on weighted samples")
  }
  if(findLambda) lambda <- NA_real_
  else {
    if(any(negative <- lambda < 0)) {
      lambda[negative] <- 0
      warning("negative value for 'lambda'; using no penalization")
    }
    lambda <- sort.int(unique(lambda), decreasing=TRUE)
  }
  if(length(lambda) == 1) crit <- "none"
  else crit <- "BIC"  # CV not yet implemented
  standardize <- isTRUE(standardize)
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
  ncores <- rep(ncores, length.out=1)
  if(is.na(ncores)) ncores <- detectCores()  # use all available cores
  if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
    ncores <- 1  # use default value
    warning("invalid value of 'ncores'; using default value")
  } else ncores <- as.integer(ncores)

  ## obtain random subsets with only 3 observations
  if(!is.null(seed)) set.seed(seed)  # set seed of random number generator
  subsets <- randomSubsets(n, 3, nsamp[1])

  ## robustly standardize data
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
  if(length(lambda) == 1) {
    fit <- fastSparseS(x=xs, y=z, findLambda=findLambda, fixedLambda=lambda,
                       initial=subsets, nistep=nistep, nkeep=nsamp[2],
                       tuning.chi=tuning.chi, bb=bb, nfpi=nfpi, tol=tol,
                       eps=eps, use.Gram=use.Gram, ncores=ncores)
  } else {
    # grid of lambda values is supplied
    names(lambda) <- seq_along(lambda)
    fit <- lapply(lambda, fastSparseS, x=xs, y=z, findLambda=FALSE,
                  initial=subsets, nistep=nistep, nkeep=nsamp[2],
                  tuning.chi=tuning.chi, bb=bb, nfpi=nfpi, tol=tol,
                  eps=eps, use.Gram=use.Gram, ncores=ncores, drop=FALSE)
    names(fit) <- names(lambda)
  }

  ## construct return object
  x <- addIntercept(x)
  if(length(lambda) == 1) {
    # backtransform coefficients, residuals and scale
    coefficients <- backtransform(fit$coefficients, centerY=muY, scaleY=sigmaY,
                                  centerX=muX, scaleX=sigmaX)
    residuals <- fit$residuals * sigmaY
    df <- modelDf(coefficients, tol)
    fit <- list(lambda=fit$lambda,
                rweights=copyNames(from=y, to=fit$weights),
                objective=fit$objective,
                coefficients=copyNames(from=x, to=coefficients),
                fitted.values=copyNames(from=y, to=y-residuals),
                residuals=copyNames(from=y, to=residuals),
                scale=fit$scale*sigmaY, df=df, standardize=standardize)
  } else {
    # backtransform coefficients, residuals and scale
    coefficients <- sapply(fit, function(f) {
      backtransform(f$coefficients, centerY=muY, scaleY=sigmaY,
                    centerX=muX, scaleX=sigmaX)
    })
    residuals <- sapply(fit, "[[", "residuals") * sigmaY
    df <- apply(coefficients, 2, modelDf, tol)
    fit <- list(lambda=lambda,
                rweights=copyNames(from=y, to=sapply(fit, "[[", "weights")),
                objective=sapply(fit, "[[", "objective"),
                coefficients=copyNames(from=x, to=coefficients),
                fitted.values=copyNames(from=y, to=y-residuals),
                residuals=copyNames(from=y, to=residuals),
                scale=sapply(fit, "[[", "scale") * sigmaY, df=df,
                standardize=standardize)
  }

  ## assign class
  class(fit) <- c("sparseS", "penModel")

  ## add information on the optimal model
  if(crit == "BIC") fit$crit <- bicSelect(fit)

  ## add center and scale estimates
  fit[c("muX", "sigmaX", "muY", "sigmaY")] <- list(muX, sigmaX, muY, sigmaY)
  ## add model data if requested
  if(isTRUE(model)) fit[c("x", "y")] <- list(x=x, y=y)

  ## return results
  fit$call <- matchedCall
  fit
}


## internal function to compute sparse S-estimator
fastSparseS <- function(fixedLambda, x, y, findLambda = FALSE, initial,
                        nistep = 10, nkeep = 10, tuning.chi = 1.547645,
                        bb = 0.5, nfpi = c(10, 500),
                        tol = .Machine$double.eps^0.5,
                        eps = .Machine$double.eps, use.Gram = TRUE,
                        ncores = 1, drop = TRUE) {
  # call C++ function
  fit <- .Call("R_fastSparseS", R_x=x, R_y=y, R_findLambda=findLambda,
               R_fixedLambda=fixedLambda, R_initial=initial, R_nistep=nistep,
               R_nkeep=nkeep, R_k=tuning.chi, R_b=bb, R_nfpi=nfpi[1],
               R_nfpiMax=nfpi[2], R_tol=tol, R_eps=eps, R_useGram=use.Gram,
               R_ncores=ncores, PACKAGE="robustHD")
  if(drop) {
    # drop the dimension of selected components
    which <- c("weights", "residuals")
    fit[which] <- lapply(fit[which], drop)
  }
  fit
}


## internal functions for expected value of M-scale under normal model
## originally written by Christophe Croux

# breakdown point belonging to constant c:   Tbsb(c, 1) / c
# value of 'bb' corresponding to constant c: Tbsb(c, 1) / c * Mrho(c, c)

# compute int(-c,c) x^(2s)dPhi(x) with Phi the cdf of N(0,1) when p=1
ksiint <- function(c, s, p) {
  (2^s) * gamma(s+p/2) * pgamma(c^2/2, s+p/2) / gamma(p/2)
}

# compute 3b/c or equivalently c*alpha (again p=1)
Tbsb <- function(c, p) {
  y1 <- ksiint(c, 1, p)*3/c - ksiint(c, 2, p)*3/(c^3) + ksiint(c, 3, p)/(c^5)
  y2 <- c*(1-pchisq(c^2, p))
  y1+y2
}

# computes c for a fixed breakdown point (p=1)
Tbsc <- function(alpha, p, maxIteration = 1000, tol = .Machine$double.eps^0.5) {
  # constant for Tukey biweight S
  c <- sqrt(qchisq(1-alpha, p))
  diff <- Inf
  i <- 1
  while ((diff > tol) && (i < maxIteration)) {
    previous <- c
    c <- Tbsb(previous, p) / alpha
    diff <- abs(previous - c)
    i <- i+1
  }
  c
}
