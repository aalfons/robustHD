# -------------------------------------------------------------------
# Author: Andreas Alfons
#         KU Leuven
#
# based on code by Jafar A. Khan, Stefan Van Aelst and Ruben H. Zamar
# -------------------------------------------------------------------

#' Robust least angle regression
#' 
#' Robustly sequence candidate predictors according to their predictive content 
#' and find the optimal model along the sequence.
#' 
#' @aliases print.rlars
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{rlars} is called.
#' @param x  a matrix or data frame containing the candidate predictors.
#' @param y  a numeric vector containing the response.
#' @param sMax  an integer vector of length two.  If a single integer is 
#' supplied, it is recycled.  The first element gives the number of predictors 
#' to be sequenced.  If it is \code{NA} (the default), predictors are sequenced 
#' as long as there are no singularity issues.  The second element gives the 
#' maximum number of predictors to be included in the final model.  If it is 
#' \code{NA} (the default), predictors may be added to the model as long as 
#' there are twice as many observations as predictors.
#' @param centerFun  a function to compute a robust estimate for the center 
#' (defaults to \code{\link[stats]{median}}).
#' @param scaleFun  a function to compute a robust estimate for the scale 
#' (defaults to \code{\link[stats]{mad}}).
#' @param const numeric; tuning constant to be used in the initial corralation 
#' estimates based on adjusted univariate winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in bivariate 
#' winsorization (defaults to 0.95).
#' @param fit  a logical indicating whether to fit submodels along the sequence 
#' (\code{TRUE}, the default) or to simply return the sequence (\code{FALSE}).
#' @param regFun  a function to compute robust linear regressions along the 
#' sequence (defaults to \code{\link[robustbase]{lmrob}}).
#' @param regArgs  a list of arguments to be passed to \code{regFun}.
#' @param crit  a character string specifying the optimality criterion to be 
#' used for selecting the final model.  Currently, only \code{"BIC"} for the 
#' Bayes information criterion is implemented.
#' @param ncores  a positive integer giving the number of processor cores to be 
#' used for parallel computing (the default is 1 for no parallelization).  If 
#' this is set to \code{NA}, all available processor cores are used.
#' @param model  a logical indicating whether the model data should be included 
#' in the returned object.
#' @param tol  a small positive numeric value.  This is used in bivariate 
#' winsorization to determine whether the initial estimate from adjusted 
#' univariate winsorization is close to 1 in absolute value.  In this case, 
#' bivariate winsorization would fail since the points form almost a straight 
#' line, and the initial estimate is returned.
#' @param \dots  additional arguments to be passed down.  For the default 
#' method, additional arguments to be passed down to 
#' \code{\link[=standardize]{robStandardize}}.
#' 
#' @return 
#' If \code{fit} is \code{FALSE}, an integer vector containing the indices of 
#' the sequenced predictors.
#'  
#' Otherwise an object of class \code{"rlars"} (inheriting from class 
#' \code{"seqModel"}) with the following components:
#' @returnItem active  an integer vector containing the indices of the 
#' sequenced predictors.
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
#' predictors.
#' @returnItem sigmaX  a numeric vector containing the standard deviations of 
#' the cleaned predictors.
#' @returnItem x  the matrix of candidate predictors (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#' 
#' @note Parallel computing for some of of the more computer-intensive 
#' computations for the inactive predictors is implemented on the C++ level 
#' via OpenMP (\url{http://openmp.org/}).
#' 
#' @author Andreas Alfons, based on code by Jafar A. Khan, Stefan Van Aelst and 
#' Ruben H. Zamar
#' 
#' @references 
#' Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust linear model 
#' selection based on least angle regression. \emph{Journal of the American 
#' Statistical Association}, \bold{102}(480), 1289--1299.
#' 
#' @seealso \code{\link{coef.seqModel}}, \code{\link{fitted.seqModel}}, 
#' \code{\link{residuals.seqModel}}, \code{\link{predict.seqModel}}, 
#' \code{\link{plot.seqModel}}
#' 
#' @example inst/doc/examples/example-rlars.R
#' 
#' @keywords regression robust
#' 
#' @export

rlars <- function(x, ...) UseMethod("rlars")


#' @rdname rlars
#' @method rlars formula
#' @export

rlars.formula <- function(formula, data, ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rlars")
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
    out <- rlars.default(x, y, ...)
    out$call <- call  # add call to return object
    out$terms <- mt   # add model terms to return object
    out
}


#' @rdname rlars
#' @method rlars default
#' @export

rlars.default <- function(x, y, sMax = NA, centerFun = median, scaleFun = mad, 
        const = 2, prob = 0.95, fit = TRUE, regFun = lmrob, regArgs = list(), 
        crit = "BIC", ncores = 1, model = TRUE, tol = .Machine$double.eps^0.5, 
        ...) {
    ## initializations
    call <- match.call()  # get function call
    call[[1]] <- as.name("rlars")
    n <- length(y)
    x <- addColnames(as.matrix(x))
    p <- ncol(x)
    sMax <- checkSMax(sMax, n, p)
    # robustly standardize data
    z <- robStandardize(y, centerFun, scaleFun, ...)   # standardize response
    xs <- robStandardize(x, centerFun, scaleFun, ...)  # standardize predictors
    # check regression function
    regControl <- getRegControl(regFun)
    regFun <- regControl$fun  # if possible, do not use formula interface
    crit <- match.arg(crit)
    ncores <- rep(ncores, length.out=1)
    if(is.na(ncores)) {
        ncores <- 0  # use all available cores
    } else if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
        ncores <- 1  # use default value
        warning("invalid value of 'ncores'; using default value")
    } else ncores <- as.integer(ncores)
    
    ## call C++ function
    active <- .Call("R_fastRlars", R_x=xs, R_y=z, R_sMax=as.integer(sMax[1]), 
        R_c=as.numeric(const), R_prob=as.numeric(prob), R_tol=as.numeric(tol), 
        scaleFun=scaleFun, R_ncores=ncores, PACKAGE="robustHD") + 1
    
    ## choose optimal model according to specified criterion
    if(isTRUE(fit)) {
        # add ones to matrix of predictors to account for intercept
        x <- addIntercept(x)
        # call function to fit models along the sequence
        s <- if(is.na(sMax[2])) NULL else 0:sMax[2]
        out <- fitModels(x, y, s=s, robust=TRUE, regFun=regFun, 
            useFormula=regControl$useFormula, regArgs=regArgs, 
            active=active, crit=crit, class="rlars")
        # add center and scale estimates
        out$muY <- attr(z, "center")
        out$sigmaY <- attr(z, "scale")
        out$muX <- attr(xs, "center")
        out$sigmaX <- attr(xs, "scale")
        if(isTRUE(model)) {
            # add model data to result
            out$x <- x
            out$y <- y
        }
        out$call <- call  # add call to return object
        out
    } else active
}
