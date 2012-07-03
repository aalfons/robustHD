# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## supplement the coefficients in a model with other useful information
## returns a data frame suitable for plotting with ggplot2

coefify <- function(model, ...) UseMethod("coefify")

coefify.seqModel <- function(model, zeros = FALSE, labels, ...) {
    # prepare coefficients and labels
    coef <- removeIntercept(t(coef(model, s=NULL)))
    sigmaX <- model$sigmaX
#    if(missing(labels)) labels <- defaultLabels(model)  # default labels
    if(!isTRUE(zeros)) {
        keep <- apply(coef != 0, 2, any)
        coef <- coef[, keep, drop=FALSE]
        sigmaX <- sigmaX[keep]
        if(!is.null(labels)) labels <- labels[keep]
    }
    # standardize coefficients
    coef <- sweep(coef, 2, sigmaX, "/", check.margin=FALSE)
    # prepare other information
    m <- ncol(coef)           # number of variables
    steps <- getSteps(model)  # step numbers
    nsteps <- length(steps)   # number of steps
    df <- model$df            # degrees of freedom
    vn <- colnames(coef)      # variable names
    # build data frame
    coefData <- data.frame(step=rep.int(steps, m), 
        df=rep.int(df, m), coefficient=as.numeric(coef), 
        variable=rep(factor(vn, levels=vn), each=nsteps))
    if(!is.null(labels)) 
        coefData$label <- rep(as.character(labels), each=nsteps)
    coefData
}

coefify.sparseLTSGrid <- function(model, fit = c("reweighted", "raw", "both"), 
        zeros = FALSE, labels, ...) {
    # prepare coefficients and labels
    fit <- match.arg(fit)
    coef <- removeIntercept(t(coef(model, s=NULL, fit=fit)))
    zeros <- isTRUE(zeros)
#    if(missing(labels)) labels <- defaultLabels(model)  # default labels
    if(!zeros) {
        keep <- apply(coef != 0, 2, any)
        coef <- coef[, keep, drop=FALSE]
        if(!is.null(labels)) labels <- labels[keep]
    }
    # check if predictor data is available to compute them
    if(is.null(x <- model$x)) {
        x <- try(model.matrix(model$terms), silent=TRUE)
        if(inherits(x, "try-error")) 
            stop("scale estimates of predictor variables not available")
    }
    x <- removeIntercept(x)
    if(!zeros) x <- x[, keep, drop=FALSE]
    # obtain scale estimates for predictors
    n <- nrow(x)
    sigmaX <- apply(x, 2, 
        function(x) {
            # standardize data
            xs <- robStandardize(x, fallback=TRUE)
            # detect good data points
            ok <- which(abs(xs) < qnorm(0.9875))
            nOk <- length(ok)
            # compute consistency factor
            if(nOk < n) {
                qn <- qnorm((nOk+n)/ (2*n))  # quantile for consistency factor
                cdelta <- 1 / sqrt(1-(2*n)/(nOk/qn)*dnorm(qn))
            } else cdelta <- 1  # consistency factor not necessary
            # compute standard deviation of good data points and multiply with 
            # consistency factor
            sd(x[ok]) * cdelta
        })
    # standardize coeffients
    coef <- sweep(coef, 2, sigmaX, "/", check.margin=FALSE)
    # prepare other information
    m <- ncol(coef)           # number of variables
    lambda <- model$lambda    # tuning parameters
    steps <- getSteps(model)  # step numbers
    sMax <- length(steps)     # number of steps
    df <- model$df            # degrees of freedom
    vn <- colnames(coef)      # variable names
    # build data frame
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        coefData <- data.frame(
            fit=rep.int(rep(factor(fits, levels=fits), each=sMax), m), 
            lambda=rep.int(lambda, 2*m), step=rep.int(steps, 2*m), 
            df=rep.int(df, 2*m), coefficient=as.numeric(coef), 
            variable=rep(factor(vn, levels=vn), each=2*sMax))
        if(!is.null(labels)) 
            coefData$label <- rep(as.character(labels), each=2*sMax)
    } else {
        coefData <- data.frame(
            lambda=rep.int(lambda, m), step=rep.int(steps, m), 
            df=rep.int(df, m), coefficient=as.numeric(coef), 
            variable=rep(factor(vn, levels=vn), each=sMax))
        if(!is.null(labels)) 
            coefData$label <- rep(as.character(labels), each=sMax)
    }
    coefData
}


#' Plot a sequence of regression models
#' 
#' Produce a plot of the coefficients or values of the optimality criterion for 
#' a sequence of regression models, such as submodels along a robust least 
#' angle regression sequence, or sparse least trimmed squares regression models 
#' for a grid of values for the penalty parameter.
#' 
#' @method plot seqModel
#' @aliases plot.rlars plot.sparseLTSGrid
#' 
#' @param x  the model fit to be plotted.
#' @param method  a character string specifying the type of plot.  Possible 
#' values are \code{"coefficients"} to plot the coefficients from the submodels 
#' via \code{\link{coefPlot}}, or \code{"crit"} to plot the values of the 
#' optimality criterion for the submodels via \code{\link{critPlot}}.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return  
#' An object of class \code{"trellis"} (see \code{\link[lattice]{xyplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{coefPlot}}, \code{\link{critPlot}}, \code{\link{rlars}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-plot.rlars.R
#' 
#' @keywords hplot
#' 
#' @export

plot.seqModel <- function(x, method = c("coefficients", "crit"), ...) {
    ## initializations
    method <- match.arg(method)
    ## call plot function
    if(method == "coefficients") coefPlot(x, ...)
    else critPlot(x, ...)
}

# ----------------------

#' Coefficient plot of a sequence of regression models
#' 
#' Produce a plot of the coefficients from a sequence of regression models, 
#' such as submodels along a robust least angle regression sequence, or sparse 
#' least trimmed squares regression models for a grid of values for the penalty 
#' parameter.
#' 
#' @aliases coefPlot.rlars
#' 
#' @param x  the model fit to be plotted.
#' @param fit  a character string specifying for which estimator to produce the 
#' plot.  Possible values are \code{"reweighted"} (the default) for the 
#' reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for both 
#' estimators.
#' @param abscissa  a character string specifying what to plot on the 
#' \eqn{x}-axis.  Possible values are \code{"step"} for the step number (the 
#' default), or \code{"df"} for the degrees of freedom.
#' @param zeros  a logical indicating whether predictors that never enter the 
#' model and thus have zero coefficients should be included in the plot 
#' (\code{TRUE}) or omitted (\code{FALSE}, the default).  This is useful if the 
#' number of predictors is much larger than the number of observations, in 
#' which case many coefficients are never nonzero.
#' @param grid  a logical indicating whether vertical grid lines should be 
#' drawn at each step.
#' @param labels  an optional character vector containing labels for the 
#' predictors.
#' @param pos  an integer position specifier for the labels.  Possible values 
#' are 1, 2, 3 and 4, respectively indicating positions below, to the left of, 
#' above and to the right of the corresponding coefficient values from the last 
#' step.
#' @param offset   a numeric value giving the offset of the labels from the 
#' corresponding coefficient values from the last step (in fractions of a 
#' character width).
#' @param \dots  for the generic function, additional arguments to be passed 
#' down to methods.  For the \code{"seqModel"} and \code{"sparseLTSGrid"} 
#' methods, additional arguments to be passed down to 
#' \code{\link[lattice]{xyplot}}.
#' 
#' @return  
#' An object of class \code{"trellis"} (see \code{\link[lattice]{xyplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[lattice]{xyplot}}, \code{\link{rlars}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-coefPlot.rlars.R
#' 
#' @keywords hplot
#' 
#' @export

coefPlot <- function(x, ...) UseMethod("coefPlot")


#' @rdname coefPlot
#' @method coefPlot seqModel
#' @export

coefPlot.seqModel <- function(x, abscissa = c("step", "df"), zeros = FALSE, 
        size = c(0.5, 2, 4), labels, offset = 1, ...) {
    ## initializations
    if(missing(labels)) labels <- defaultLabels(x)  # default labels
    ## extract coefficient data extended with other information
    coefData <- coefify(x, zeros=zeros, labels=labels)
    ## construct data frame for labels
    maxStep <- max(coefData$step)
    labelData <- coefData[coefData$step == maxStep, ]
    ## call workhorse function
    ggCoefPlot(coefData, labelData, abscissa=abscissa, size=size, 
        offset=offset, ...)
}


#' @rdname coefPlot
#' @method coefPlot sparseLTSGrid
#' @export

coefPlot.sparseLTSGrid <- function(x, fit = c("reweighted", "raw", "both"), 
        abscissa = c("step", "df"), zeros = FALSE, size = c(0.5, 2, 4), 
        labels, offset = 1, ...) {
    ## initializations
    fit <- match.arg(fit)
    abscissa <- match.arg(abscissa)
    if(missing(labels)) labels <- defaultLabels(x)  # default labels
    ## extract coefficient data extended with other information
    coefData <- coefify(x, fit=fit, zeros=zeros, labels=labels)
    ## construct data frame for labels
    maxX <- max(coefData[, abscissa])
    labelData <- coefData[coefData[, abscissa] == maxX, ]
    if(abscissa == "df") {
        # maximum degree of freedom may occur in more than one step
        # ensure that label is only drawn once for largest step number
        by <- if(fit == "both") c("fit", "variable") else "variable"
        keep <- split(rownames(labelData), labelData[, by])
        keep <- sapply(keep, tail, 1)
        labelData <- labelData[keep, ]
    }
    ## call workhorse function
    p <- ggCoefPlot(coefData, labelData, abscissa=abscissa, size=size, 
        offset=offset, ...)
    if(fit == "both") {
        f <- as.formula(paste(".", "fit", sep="~"))
        p <- p + facet_grid(f)
    }
    p
}


## workhorse function
ggCoefPlot <- function(coefData, labelData, abscissa = c("step", "df"), 
        zeros = FALSE, size = c(0.5, 2, 4), labels, offset = 1, ..., 
        mapping, data, xlab, ylab) {
    # initializations
    abscissa <- match.arg(abscissa)
    size <- as.numeric(size)
    size <- c(size, rep.int(NA, max(0, 3-length(size))))[1:3]  # ensure length 3
    size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
    # define default axis labels
    if(missing(xlab)) 
        xlab <- switch(abscissa, step="Step", df="Degrees of freedom")
    if(missing(ylab)) ylab <- "Standardized coefficients"
    # define aesthetic mapping for plotting coefficients
    coefMapping <- aes_string(x=abscissa, y="coefficient", color="variable")
    # define aesthetic mapping for plotting x-axis grid and labels
    offset <- paste(rep.int(" ", offset), collapse="")  # whitespace
    labelData$label <- paste(offset, labelData$label, sep="")
    labelMapping <- aes_string(x=abscissa, y="coefficient", label="label")
    # draw minor grid lines for each step, but leave major grid lines and tick 
    # marks pretty
    gridX <- unique(coefData[, abscissa])
    # create plot
    ggplot(coefData) + 
        geom_line(coefMapping, size=size[1], ...) + 
        geom_point(coefMapping, size=size[2], ...) + 
        geom_text(labelMapping, data=labelData, 
            hjust=0, size=size[3], alpha=0.3) + 
        scale_x_continuous(minor_breaks=gridX) + 
        opts(legend.position="none") + 
        labs(x=xlab, y=ylab)
}

# ----------------------

#' Optimality criterion plot of a sequence of regression models
#' 
#' Produce a plot of the values of the optimality criterion for a sequence of 
#' regression models, such as submodels along a robust least angle regression 
#' sequence, or sparse least trimmed squares regression models for a grid of 
#' values for the penalty parameter.
#' 
#' @aliases critPlot.rlars
#' 
#' @param x  the model fit to be plotted.
#' @param fit  a character string specifying for which estimator to produce the 
#' plot.  Possible values are \code{"reweighted"} (the default) for the 
#' reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for both 
#' estimators.
#' @param \dots  for the generic function, additional arguments to be passed 
#' down to methods.  For the \code{"seqModel"} and \code{"sparseLTSGrid"} 
#' methods, additional arguments to be passed down to 
#' \code{\link[lattice]{xyplot}}.
#' 
#' @return  
#' An object of class \code{"trellis"} (see \code{\link[lattice]{xyplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[lattice]{xyplot}}, \code{\link{rlars}}, 
#' \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-critPlot.rlars.R
#' 
#' @keywords hplot
#' 
#' @export

critPlot <- function(x, ...) UseMethod("critPlot")


#' @rdname critPlot
#' @method critPlot seqModel
#' @export

critPlot.seqModel <- function(x, ...) {
    ## extract information from object
    crit <- x$crit
    ## build data.frame for lattice graphics
    critData <- data.frame(getSteps(x), x$critValues)
    names(critData) <- c("Step", crit)
    ## define local version of 'xyplot' with different default values
    ## this also avoids error message if 'data' argument is supplied
    localXyplot <- function(..., data, type) {
        if(missing(type)) type <- "b"
        xyplot(..., data=critData, type=type)
    }
    ## call 'xyplot
    form <- as.formula(paste(crit, "Step", sep=" ~ "))  # formula
    localXyplot(form, ...)
}


#' @rdname critPlot
#' @method critPlot sparseLTSGrid
#' @export

critPlot.sparseLTSGrid <- function(x, fit = c("reweighted", "raw", "both"), 
        ...) {
    ## initializations
    fit <- match.arg(fit)
    ## extract information from object
    lambda <- x$lambda
    crit <- x$crit
    ## build data.frame for lattice graphics
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        sMax <- length(lambda)
        critData <- data.frame(factor(rep(fits, each=sMax), levels=fits), 
            rep.int(lambda, 2), c(x$critValues, x$raw.critValues))
        names(critData) <- c("Fit", "lambda", crit)
    } else {
        critValues <- if(fit == "reweighted") x$critValues else x$raw.critValues
        critData <- data.frame(lambda, critValues)
        names(critData) <- c("lambda", crit)
    }
    ## define local version of 'xyplot' with different default values
    ## this also avoids error message if 'data' argument is supplied
    localXyplot <- function(..., data, type) {
        if(missing(type)) type <- "b"
        xyplot(..., data=critData, type=type)
    }
    ## call 'xyplot
    conditional <- if(fit == "both") "Fit" else NULL
    form <- getFormula(crit, "lambda", conditional)  # formula
    localXyplot(form, ...)
}

# ----------------------

#' Diagnostic plots for sparse LTS regression models
#' 
#' Produce diagnostic plots for sparse least trimmed squares regression 
#' models.  Four plots are currently implemented.
#' 
#' In the normal Q-Q plot of the standardized residuals, a reference line is 
#' drawn through the first and third quartile.  The \code{id.n} observations 
#' with the largest distances from that line are identified by a label (the 
#' observation number).  The default for \code{id.n} is the number of 
#' regression outliers, i.e., the number of observations with outlier weight 
#' equal to 0  (see \code{\link[=weights.sparseLTS]{weights}}).
#' 
#' In the plots of the standardized residuals versus their index or the fitted 
#' values, horizontal reference lines are drawn at 0 and +/-2.5.  The 
#' \code{id.n} observations with the largest absolute values of the 
#' standardized residuals are identified by a label (the observation 
#' number).  The default for \code{id.n} is the number of regression outliers, 
#' i.e., the number of observations with outlier weight equal to 0  (see 
#' \code{\link[=weights.sparseLTS]{weights}}).
#' 
#' For the regression diagnostic plot, the robust Mahalanobis distances of the 
#' predictor variables are computed via the MCD based on only those predictors 
#' with non-zero coefficients.  Horizontal reference lines are drawn at +/-2.5 
#' and a vertical reference line is drawn at the upper 97.5\% quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution with \eqn{p} degrees of freedom, 
#' where \eqn{p} denotes the number of predictors with non-zero 
#' coefficients.  The \code{id.n} observations with the largest absolute values 
#' of the standardized residuals and/or largest robust Mahalanobis distances 
#' are identified by a label (the observation number).  The default for 
#' \code{id.n} is the number of all outliers: regression outliers (i.e., 
#' observations with outlier weight equal to 0, see 
#' \code{\link[=weights.sparseLTS]{weights}}) and leverage points (i.e., 
#' observations with robust Mahalanobis distance larger than the 97.5\% 
#' quantile of the \eqn{\chi^{2}}{chi-squared} distribution with \eqn{p} 
#' degrees of freedom).
#' 
#' @param x  the model fit for which to produce diagnostic plots.
#' @param fit  a character string specifying for which fit to produce 
#' diagnostic plots.  Possible values are \code{"reweighted"} (the default) for 
#' diagnostic plots for the reweighted fit, \code{"raw"} for diagnostic plots 
#' for the raw fit, or \code{"both"} for diagnostic plots for both fits.
#' @param which  a character string indicating which plot to show.  Possible 
#' values are \code{"all"} (the default) for all of the following, \code{"rqq"} 
#' for a normal Q-Q plot of the standardized residuals, \code{"rindex"} for a 
#' plot of the standardized residuals versus their index, \code{"rfit"} for a 
#' plot of the standardized residuals versus the fitted values, or 
#' \code{"rdiag"} for a regression diagnostic plot  (standardized residuals 
#' versus robust Mahalanobis distances of the predictor variables).
#' @param ask  a logical indicating whether the user should be asked before 
#' each plot (see \code{\link[grDevices]{devAskNewPage}}). The default is to 
#' ask if all plots are requested and not ask otherwise.
#' @param id.n  an integer giving the number of the most extreme observations 
#' to be identified by a label.  The default is to use the number of identified 
#' outliers, which can be different for the different plots.  See 
#' \dQuote{Details} for more information.
#' @param \dots  for the generic function \code{diagnosticPlot}, additional 
#' arguments to be passed down to methods.  For the \code{"sparseLTSGrid"} 
#' method of \code{diagnosticPlot}, additional arguments to be passed down to 
#' the \code{"sparseLTS"} method.  For the \code{"sparseLTS"} method of 
#' \code{diagnosticPlot}, additional arguments to be passed down to 
#' \code{\link[lattice]{xyplot}}.  For the \code{"sparseLTS"} method of 
#' \code{plot}, additional arguments to be passed down to \code{diagnosticPlot}.
#' 
#' @return  
#' If only one plot is requested, an object of class \code{"trellis"} (see 
#' \code{\link[lattice]{xyplot}}), otherwise a list of such objects.
#' 
#' @author Andreas Alfons, partly based on code by Valentin Todorov
#' 
#' @seealso \code{\link[graphics]{plot}}, 
#' \code{\link[robustbase:ltsPlot]{plot.lts}}, 
#' \code{\link{sparseLTS}}, \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-diagnosticPlot.R
#' 
#' @keywords hplot
#' 
#' @export
#' @import robustbase

diagnosticPlot <- function(x, ...) UseMethod("diagnosticPlot")


#' @rdname diagnosticPlot
#' @method diagnosticPlot sparseLTS
#' @export

diagnosticPlot.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"), 
        which = c("all", "rqq","rindex", "rfit", "rdiag"),
        ask = (which == "all"), id.n = NULL, ...) {
    # initializations
    fit <- match.arg(fit)
    which <- match.arg(which)
    weights <- as.matrix(weights(x, fit=fit))
    if(id.n.default <- is.null(id.n)) {
        if(fit != "rdiag") {
            id.n <- apply(weights, 2, function(w) length(which(w == 0)))
        }
    } else {
        d <- dim(weights)
        id.n <- rep(as.integer(id.n), length.out=d[2])
        if(!isTRUE(all(id.n >= 0 & id.n <= d[1]))) {
            stop("'id.n' must be in {1,..,", d[1], "}")
        }
    }
    scale <- switch(fit, reweighted=x$scale, raw=x$raw.scale, 
        both=c(reweighted=x$scale, raw=x$raw.scale))
    if(all(scale <= 0)) {
        stop("plots not available (residual scale equal to 0)")
    }
    # call functions for selected plots
    if(which == "all") {
        devAskNewPage(ask)  # ask for new page (if requested)
        # Q-Q plot
        tmp <- try(qqmath(x, fit=fit, id.n=id.n, ...), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res <- list(rqq=tmp)
        }
        # residuals vs index plot
        tmp <- try(indexplot(x, fit=fit, id.n=id.n, ...), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res$rindex <- tmp
        }
        # residuals vs fitted plot
        tmp <- try(fitplot(x, fit=fit, id.n=id.n, ...), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res$rfit <- tmp
        }
        # regression diagnostic plot
        if(id.n.default) id.n <- NULL
        tmp <- try(diagplot(x, fit=fit, id.n=id.n, ...), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res$rdiag <- tmp
        }
    } else if(which == "rqq") {
        # Q-Q plot
        res <- qqmath(x, fit=fit, id.n=id.n, ...)
        print(res)
    } else if(which == "rindex") {
        # residuals vs index plot
        res <- indexplot(x, fit=fit, id.n=id.n, ...)
        print(res)
    } else if(which == "rfit") {
        # residuals vs fitted plot
        res <- fitplot(x, fit=fit, id.n=id.n, ...)
        print(res)
    } else if(which == "rdiag") {
        # regression diagnostic plot
        res <- diagplot(x, fit=fit, id.n=id.n, ...)
        print(res)
    }
    invisible(res)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot sparseLTSGrid
#' @export

# TODO: allow to select which submodels to plot
diagnosticPlot.sparseLTSGrid <- function(x, ...) {
    # initializations
    s <- x$sOpt
    raw.s <- x$raw.sOpt
    lambda <- unname(x$lambda)
    if(s == raw.s) {
        lambda <- lambda[s]
    } else lambda <- c(reweighted=lambda[s], raw=lambda[raw.s])
    # this is a bit of a hack: create an object of class "sparseLTS" that 
    # contains the information of the optimal model for the reweighted 
    # estimator and the optimal model for the raw estimator
    x$coefficients <- coef(x, s=s, fit="reweighted")
    x$fitted.values <- fitted(x, s=s, fit="reweighted")
    x$residuals <- residuals(x, s=s, fit="reweighted")
    x$weights <- weights(x, s=s, fit="reweighted")
    x$raw.coefficients <- coef(x, s=raw.s, fit="raw")
    x$raw.residuals <- residuals(x, s=raw.s, fit="raw")
    x$raw.weights <- weights(x, s=raw.s, fit="raw")
    x$best <- x$best[, raw.s]
    x$objective <- unname(x$objective[raw.s])
    x$center <- unname(x$center[s])
    x$scale <- unname(x$scale[s])
    x$lambda <- lambda
    x$cnp2 <- unname(x$cnp2[s])
    x$df <- unname(x$df[s])
    x$raw.center <- unname(x$raw.center[raw.s])
    x$raw.scale <- unname(x$raw.scale[raw.s])
    x$raw.cnp2 <- unname(x$raw.cnp2[raw.s])
    x$crit <- x$critValues <- x$sOpt <- x$raw.critValues <- x$raw.sOpt <- NULL
    class(x) <- "sparseLTS"
    # call the method for class "sparseLTS"
    diagnosticPlot(x, ...)
}


#' @rdname diagnosticPlot
#' @method plot sparseLTS
#' @export

plot.sparseLTS <- function(x, ...) diagnosticPlot(x, ...)


## Q-Q plots for sparse LTS residuals
qqmath.sparseLTS <- function(x, data, fit = c("reweighted", "raw", "both"), 
        id.n = NULL, ...) {
    # initializations
    fit <- match.arg(fit)
    weights <- as.matrix(weights(x, fit=fit))
    if(is.null(id.n)) {
        id.n <- apply(weights, 2, function(w) length(which(w == 0)))
    } else {
        d <- dim(weights)
        id.n <- rep(as.integer(id.n), length.out=d[2])
        if(!isTRUE(all(id.n >= 0 & id.n <= d[1]))) {
            stop("'id.n' must be in {1,..,", d[1], "}")
        }
    }    
    # construct data frame in lattice format and call internal function
    qq <- getLatticeDataRqq(x, fit=fit)
    localQqmath(qq, id.n=id.n, ...)
}

# internal function for Q-Q plots
localQqmath <- function(qq, panel = panel.rqq, id.n, 
        main = "Normal Q-Q plot", 
        xlab = "Quantiles of the standard normal distribution", 
        ylab = "Standardized sparse LTS residual", ..., 
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups, f.value, distribution, tails.n) {
    # construct formula for call to xyplot()
    conditional <- if("Fit" %in% names(qq)) "Fit" else NULL
    f <- getFormula(NULL, "Residual", conditional)
    # call qqmath()
    qqmath(f, data=qq, panel=panel, id.n=id.n, 
        main=main, xlab=xlab, ylab=ylab, ...)
}

# panel function for Q-Q plots
panel.rqq <- function(x, id.n, ..., identifier = "qqmath") {
    # draw a line through the first and third quartiles and use intercept and 
    # slope to order the observations according to their distance from the line
    l <- panel.rqqline(x, ...)
    # create Q-Q plot as in function qqnorm() from package 'stats'
    qq <- qqnorm(x, plot.it=FALSE, datax=FALSE)
    panel.xyplot(qq$x, qq$y, ..., identifier=identifier)
    # plot labels for observations with largest distance from the line
    d <- abs(qq$y - l$intercept - l$slope * qq$x)
    panel.label(qq$x, qq$y, ord=d, id.n=id.n[packet.number()], ...)
}

# panel function to add line through the first and third quartiles
panel.rqqline <- function(x, qtype = 7, ...) {
    y <- quantile(x, c(0.25, 0.75), na.rm=TRUE, names=FALSE, type=qtype)
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y) / diff(x)
    intercept <- y[1] - slope * x[1]
    if(is.finite(intercept) && is.finite(slope)) {
        panel.refline(a=intercept, b=slope)
    }
    invisible(list(intercept=intercept, slope=slope))
}


## plot sparse LTS residuals against their index
indexplot <- function(x, fit = c("reweighted", "raw", "both"), 
        panel = panel.residuals, id.n, main = "Residuals vs index", 
        xlab = "Index", ylab = "Standardized sparse LTS residual", ..., 
        # the following arguments are defined so that they aren't supplied twice
        formula, data, groups) {
    # construct data frame in lattice format
    xy <- getLatticeDataIndex(x, fit=fit)
    # construct formula for call to xyplot()
    conditional <- if("Fit" %in% names(xy)) "Fit" else NULL
    f <- getFormula("Residual", "Index", conditional)
    # call xyplot()
    xyplot(f, data=xy, panel=panel, id.n=id.n, 
        main=main, xlab=xlab, ylab=ylab, ...)
}


## plot sparse LTS residuals against fitted values
fitplot <- function(x, fit = c("reweighted", "raw", "both"), 
        panel = panel.residuals, id.n, main = "Residuals vs fitted values", 
        xlab = "Fitted values", ylab = "Standardized sparse LTS residual", ..., 
        # the following arguments are defined so that they aren't supplied twice
        formula, data, groups) {
    # construct data frame in lattice format
    xy <- getLatticeDataFit(x, fit=fit)
    # construct formula for call to xyplot()
    conditional <- if("Fit" %in% names(xy)) "Fit" else NULL
    f <- getFormula("Residual", "Fitted", conditional)
    # call xyplot()
    xyplot(f, data=xy, panel=panel, id.n=id.n, 
        main=main, xlab=xlab, ylab=ylab, ...)
}

# panel function to plot sparse LTS residuals against index or fitted values
panel.residuals <- function(x, y, id.n, ...) {
    # plot horizontal reference lines
    panel.refline(h=-2.5)
    panel.refline(h=0)
    panel.refline(h=2.5)
    # plot residuals against index or fitted values
    panel.xyplot(x, y, ...)
    # plot labels for most extreme observations
    panel.label(x, y, ord=abs(y), id.n=id.n[packet.number()], ...)
}


## regression diagnostic plot
diagplot <- function(x, fit = c("reweighted", "raw", "both"), 
        panel = panel.diag, id.n = NULL, 
        main = "Regression diagnostic plot", 
        xlab = "Robust distance computed by MCD", 
        ylab = "Standardized sparse LTS residual", ..., 
        # the following arguments are defined so that they aren't supplied twice
        formula, data, groups, p) {
    ## initializations
    object <- x
    # extract predictor matrix
    terms <- delete.response(object$terms)  # extract terms for model matrix
    if(is.null(x <- object$x)) {
        x <- try(model.matrix(terms), silent=TRUE)
        if(inherits(x, "try-error")) stop("model data not available")
    }
    if(object$intercept) x <- removeIntercept(x)
    n <- nrow(x)
    # extract coefficients
    fit <- match.arg(fit)
    coefficients <- as.matrix(coef(object, fit=fit))
    if(object$intercept) coefficients <- coefficients[-1, , drop=FALSE]
    if(fit != "both") colnames(coefficients) <- fit
    significant <- apply(coefficients, 2, function(x) x != 0)
    p <- apply(significant, 2, sum)
    if(all(p <= 0)) {
        stop("all coefficients equal to 0")
    }
    # extract residuals and outlier weights
    residuals <- as.matrix(residuals(object, fit=fit, standardized=TRUE))
    weights <- as.matrix(weights(object, fit=fit))
    ## compute MCD distances using significant variables
    # adjust alpha since MCD computes subset size depending on n and p
    h <- object$quan
    n2 <- (n+p+1) %/% 2
    alpha <- pmin((h - 2*n2 + n) / (2 * (n - n2)), 1)
    # compute MCD distances
    if(fit == "reweighted" || fit == "both") {
        # check fraction for subset size
        if(alpha[1] < 0.5) {
            alpha[1] <- 0.5
            warning(sprintf("cannot compute MCD with h = %d; using h = %d", 
                    object$quan, h.alpha.n(alpha[1], n, p[1])))
        }
        # compute distances
        RD <- try({
                xs <- x[, significant[, 1], drop=FALSE]
                mcd <- covMcd(xs, alpha=alpha[1])
                sqrt(mahalanobis(xs, mcd$center, mcd$cov))
            })
        if(inherits(RD, "try-eror")) {
            msg <- "robust distances cannot be computed"
            if(fit == "both") {
                RD <- rep.int(NA, nrow(residuals))
                warning(msg)
            } else stop(msg)
        }
        RD <- as.matrix(RD)
    } else RD <- NULL
    if(fit == "raw" || fit == "both") {
        # check fraction for subset size
        if(alpha["raw"] < 0.5) {
            alpha["raw"] <- 0.5
            warning(sprintf("cannot compute MCD with h = %d; using h = %d", 
                    object$quan, h.alpha.n(alpha["raw"], n, p["raw"])))
        }
        # compute distances
        tmp <- try({
                xs <- x[, significant[, "raw"], drop=FALSE]
                mcd <- covMcd(xs, alpha=alpha["raw"])
                sqrt(mahalanobis(xs, mcd$center, mcd$cov))
            })
        if(inherits(tmp, "try-eror")) {
            msg <- "robust distances cannot be computed"
            if(fit == "both") {
                tmp <- rep.int(NA, nrow(residuals))
                warning(msg)
            } else stop(msg)
        }
        RD <- cbind(RD, tmp)
    }
    colnames(RD) <- colnames(residuals)
    rownames(RD) <- rownames(x)
    # construct data frame in lattice format
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        xy <- data.frame(Fit=factor(rep(fits, each=nrow(residuals)), levels=fits), 
            RD=c(RD[, "reweighted"], RD[, "raw"]),
            Residual=c(residuals[, "reweighted"], residuals[, "raw"]))
    } else {
        xy <- data.frame(RD=RD, Residual=residuals)
    }
    # construct formula for call to xyplot()
    conditional <- if("Fit" %in% names(xy)) "Fit" else NULL
    f <- getFormula("Residual", "RD", conditional)
    # call xyplot()
    xyplot(f, data=xy, panel=panel, weights=weights, df=p, 
        id.n=id.n, main=main, xlab=xlab, ylab=ylab, ...)
}

# panel function for regression diagnostic plot
panel.diag <- function(x, y, weights, df, id.n = NULL, ...) {
    # plot horizontal reference lines
    panel.refline(h=-2.5)
    panel.refline(h=2.5)
    # plot vertical reference lines
    i <- packet.number()
    q <- sqrt(qchisq(0.975, df[i]))
    panel.refline(v=max(q, 2.5))
    # plot standardized residuals against robust distances
    panel.xyplot(x, y, ...)
    # plot labels for most extreme observations
    ord <- pmax.int(abs(x/2.5), abs(y/q))
    if(is.null(id.n)) {
        id.n <- length(unique(c(which(x > q), which(weights[, i] == 0))))
    } else id.n <- id.n[i]
    panel.label(x, y, ord=ord, id.n=id.n, ...)
}


# ----------------------

## utilities for plot functions

# get formula for plot functions
getFormula <- function(left, right, conditional = NULL) {
    if(is.null(conditional)) {
        as.formula(paste(left, "~", right))
    } else as.formula(paste(left, "~", right, "|", conditional))
}

# get data in the correct format for Q-Q plots
getLatticeDataRqq <- function(x, fit = c("reweighted", "raw", "both")) {
    residuals <- residuals(x, fit=fit, standardized=TRUE)
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        data.frame(Fit=factor(rep(fits, each=nrow(residuals)), levels=fits), 
            Residual=c(residuals[, "reweighted"], residuals[, "raw"]))
    } else data.frame(Residual=residuals)
}

# get data in the correct format for plotting residuals against their index
getLatticeDataIndex <- function(x, fit = c("reweighted", "raw", "both")) {
    residuals <- residuals(x, fit=fit, standardized=TRUE)
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        data.frame(Fit=factor(rep(fits, each=nrow(residuals)), levels=fits), 
            Index=rep.int(seq_len(nrow(residuals)), 2),
            Residual=c(residuals[, "reweighted"], residuals[, "raw"]))
    } else data.frame(Index=seq_along(residuals), Residual=residuals)
}

# get data in the correct format for plotting residuals against fitted values
getLatticeDataFit <- function(x, fit = c("reweighted", "raw", "both")) {
    fitted <- try(fitted(x, fit=fit), silent=TRUE)
    residuals <- residuals(x, fit=fit, standardized=TRUE)
    if(fit == "both") {
        # check if fitted values of raw estimator are available
        if(is.null(dim(fitted))) {
            fitted <- cbind(reweighted=fitted, raw=rep.int(NA, length(fitted)))
        }
        fits <- c("reweighted", "raw")
        # construct data frame
        data.frame(Fit=factor(rep(fits, each=nrow(residuals)), levels=fits), 
            Fitted=c(fitted[, "reweighted"], fitted[, "raw"]),
            Residual=c(residuals[, "reweighted"], residuals[, "raw"]))
    } else {
        # check if fitted values of raw estimator are available
        if(inherits(fitted, "try-error")) {
            # throw warning containing the error message
            warn <- gsub("Error in", "In", fitted)
            warning(warn, call.=FALSE)
            fitted <- rep.int(NA, length(residuals))
        }
        # construct data frame
        data.frame(Fitted=fitted, Residual=residuals)
    }
}

# panel function to add labels for points with large distances
panel.label <- function(x, y, ord, lab, id.n, ...) {
    if(id.n > 0) {
        n <- length(y)
        which <- order(ord)[(n - id.n + 1):n]
        lab <- if(missing(lab)) which else lab[which]
        ## how to adjust the labels?
        ## a) pos=4 (to the left of the observation)
        ## b) additionaly to pos specify offset=0.2 (fraction of a character)
        panel.text(x[which], y[which], labels=lab, pos = 4, offset = 0.2, ...)
    }
}
