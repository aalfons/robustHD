# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

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
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
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
#' @param size  a numeric vector of length three giving the line width, the 
#' point size and the label size, respectively.
#' @param labels  an optional character vector containing labels for the 
#' predictors.  Plotting labels can be suppressed by setting this to 
#' \code{NULL}.
#' @param offset   an integer giving the offset of the labels from the 
#' corresponding coefficient values from the last step (i.e., the number of 
#' blank characters to be prepended to the label).
#' @param \dots  for the generic function, additional arguments to be passed 
#' down to methods.  For the \code{"seqModel"} and \code{"sparseLTSGrid"} 
#' methods, additional arguments to be passed down to 
#' \code{\link[ggplot2]{geom_line}} and \code{\link[ggplot2]{geom_point}}.
#' 
#' @return  
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}}, 
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
    # draw minor grid lines for each step, but leave 
    # major grid lines and tick marks pretty
    gridX <- unique(coefData[, abscissa])
    # create plot
    ggplot(coefData) + 
        geom_line(coefMapping, size=size[1], ...) + 
        geom_point(coefMapping, size=size[2], ...) + 
        geom_text(labelMapping, data=labelData, 
            hjust=0, size=size[3], alpha=0.4) + 
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
#' @param size  a numeric vector of length two giving the line width and the 
#' point size, respectively.
#' @param \dots  for the generic function, additional arguments to be passed 
#' down to methods.  For the \code{"seqModel"} and \code{"sparseLTSGrid"} 
#' methods, additional arguments to be passed down to 
#' \code{\link[ggplot2]{geom_line}} and \code{\link[ggplot2]{geom_point}}.
#' 
#' @return  
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}}, 
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

critPlot.seqModel <- function(x, size = c(0.5, 2), ...) {
    ## extract information from object
    crit <- x$crit
    ## construct data frame for ggplot2 graphics
    critData <- data.frame(getSteps(x), x$critValues)
    names(critData) <- c("step", crit)
    ## call workhorse function
    ggCritPlot(critData, abscissa="step", size=size, ...)
}


#' @rdname critPlot
#' @method critPlot sparseLTSGrid
#' @export

critPlot.sparseLTSGrid <- function(x, fit = c("reweighted", "raw", "both"), 
        size = c(0.5, 2), ...) {
    ## initializations
    fit <- match.arg(fit)
    ## extract information from object
    lambda <- x$lambda
    crit <- x$crit
    ## construct data frame for ggplot2 graphics
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        sMax <- length(lambda)
        critData <- data.frame(rep(factor(fits, levels=fits), each=sMax), 
            rep.int(lambda, 2), c(x$critValues, x$raw.critValues))
        names(critData) <- c("fit", "lambda", crit)
    } else {
        critValues <- if(fit == "reweighted") x$critValues else x$raw.critValues
        critData <- data.frame(lambda, critValues)
        names(critData) <- c("lambda", crit)
    }
    ## call workhorse function
    p <- ggCritPlot(critData, abscissa="lambda", size=size, ...)
    if(fit == "both") {
        f <- as.formula(paste(".", "fit", sep="~"))
        p <- p + facet_grid(f)
    }
    p
}


## workhorse function
ggCritPlot <- function(critData, abscissa = c("step", "lambda"), 
        size = c(0.5, 2), ..., mapping, data, xlab, ylab) {
    # initializations
    abscissa <- match.arg(abscissa)
    crit <- setdiff(names(critData), c("fit", "step", "lambda"))
    size <- as.numeric(size)
    size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
    size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
    # define default axis labels
    if(missing(xlab)) xlab <- switch(abscissa, step="Step", lambda="lambda")
    if(missing(ylab)) ylab <- crit
    # define aesthetic mapping for plotting coefficients
    mapping <- aes_string(x=abscissa, y=crit)
    # draw minor grid lines for each step, but leave 
    # major grid lines and tick marks pretty
    gridX <- unique(critData[, abscissa])
    # create plot
    ggplot(critData, mapping) + 
        geom_line(size=size[1], ...) + 
        geom_point(size=size[2], ...) + 
        scale_x_continuous(minor_breaks=gridX) + 
        labs(x=xlab, y=ylab)
}

# ----------------------

## supplement the fitted values and residuals of a model with other useful 
## information for diagnostic plots
## argument 'data' is currently ignored
## returns a data frame suitable for plotting with ggplot2

fortify.sparseLTS <- function(model, data, 
        fit = c("reweighted", "raw", "both"), ...) {
    ## initializations
    fit <- match.arg(fit)
    ## construct data frame with all information for plotting
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        ## recursive call for each fit
        reweighted <- fortify(model, fit="reweighted", ...)
        raw <- fortify(model, fit="raw", ...)
        ## combine data for Q-Q reference line
        qql <- data.frame(fit=factor(fits, levels=fits), 
            rbind(attr(reweighted, "qqLine"), attr(raw, "qqLine")), 
            row.names=NULL)
        ## combine data for cutoff chi-squared quantile
        q <- data.frame(fit=factor(fits, levels=fits), 
            rbind(attr(reweighted, "q"), attr(raw, "q")), 
            row.names=NULL)
        ## combine results
        n <- c(nrow(reweighted), nrow(raw))
        data <- data.frame(fit=rep.int(factor(fits, levels=fits), n), 
            rbind(reweighted, raw), row.names=NULL)
        attr(data, "qqLine") <- qql
        attr(data, "q") <- q
    } else {
        ## extract fitted values
        fitted <- fitted(model, fit=fit)
        ## extract standardized residuals
        residuals <- residuals(model, fit=fit, standardized=TRUE)
        n <- length(residuals)  # number of observations
        ## extract outlier weights
        weights <- weights(model, fit=fit)
        ## compute theoretical quantiles and distances from Q-Q reference line
        theoretical <- qqNorm(residuals)
        qql <- qqLine(residuals)  # Q-Q reference line
        qqd <- abs(residuals - qql$intercept - qql$slope * theoretical)
        ## compute MCD distances using significant variables
        # extract predictor matrix
        terms <- delete.response(model$terms)  # extract terms for model matrix
        ok <- TRUE
        if(is.null(x <- model$x)) {
            x <- try(model.matrix(terms), silent=TRUE)
            if(inherits(x, "try-error")) {
                ok <- FALSE
                warning("model data not available")
            }
        }
        if(ok) {
            if(model$intercept) x <- removeIntercept(x)
            # extract coefficients
            coefficients <- coef(model, fit=fit)
            if(model$intercept) coefficients <- removeIntercept(coefficients)
            significant <- which(coefficients != 0)
            p <- length(significant)
            if(p == 0) {
                ok <- FALSE
                warning("all coefficients equal to 0")
            }
        }
        if(ok) {
            # adjust alpha since MCD computes subset size depending on n and p
            h <- model$quan
            n2 <- (n+p+1) %/% 2
            alpha <- pmin((h - 2*n2 + n) / (2 * (n - n2)), 1)
            # check fraction for subset size
            if(alpha < 0.5) {
                alpha <- 0.5
                warning(sprintf("cannot compute MCD with h = %d; using h = %d", 
                        model$quan, h.alpha.n(alpha, n, p)))
            }
            # compute distances
            rd <- try({
                    xs <- x[, significant, drop=FALSE]
                    mcd <- covMcd(xs, alpha=alpha)
                    if(fit == "reweighted") {
                        center <- mcd$center
                        cov <- mcd$cov
                    } else {
                        center <- mcd$raw.center
                        cov <- mcd$raw.cov
                    }
                    sqrt(mahalanobis(xs, center, cov))
                })
            if(inherits(rd, "try-eror")) {
                ok <- FALSE
                warning("robust distances cannot be computed")
            }
        }
        if(!ok) rd <- rep.int(NA, n)
        # take maximum of the distances in the x- and y-space, divided by the 
        # respective other cutoff point
        q <- sqrt(qchisq(0.975, p))
        xyd <- pmax.int(abs(rd/2.5), abs(residuals/q))
        ## construct indicator variables for leverage points
        leverage <- rd > q
        ## classify data points
        class <- ifelse(weights == 0, "outlier", "good")
        class <- factor(class, levels=c("outlier", "good"))
        ## construct data frame
        data <- data.frame(index=seq_len(n), fitted=fitted, residual=residuals, 
            theoretical=theoretical, qqd=qqd, rd=rd, xyd=xyd, weight=weights, 
            leverage=leverage, classification=class)
        attr(data, "qqLine") <- as.data.frame(qql)
        attr(data, "q") <- data.frame(q=max(q, 2.5))
    }
    ## return data frame
    data
}


## construct data frame for labels based on some order
labelify <- function(data, which, id.n = NULL) {
    # initializations
    if(isTRUE(id.n < 1)) return(NULL)
    by <- intersect(names(data), "fit")
    ord <- data[, which]
    if(which == "residual") ord <- abs(ord)
    if(is.null(id.n)) {
        # define outlier indicator
        out <- data[, "weight"] == 0                       # regression outliers
        if(which == "xyd") out <- out | data[, "leverage"] # all outliers
    }
    # find the id.n largest observations to keep
    # if NULL, id.n is computed as the number of outliers
    if(length(by) == 0) {
        # use the whole data set
        if(is.null(id.n)) id.n <- sum(out)
        if(id.n == 0) return(NULL)
        keep <- head(order(ord, decreasing=TRUE), id.n)
    } else {
        # split the data set according to the selected variables
        keep <- tapply(seq_len(nrow(data)), data[, by, drop=FALSE], 
            function(i) {
                if(is.null(id.n)) id.n <- sum(out[i])
                if(id.n > 0) {
                    largest <- head(order(ord[i], decreasing=TRUE), id.n)
                    i[largest]
                }
            })
        # combine indices to keep
        keep <- unlist(keep, use.names=FALSE)
        if(length(keep) == 0) return(NULL)
    }
    # return data frame with selected observations
    data[keep, ]
}


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
#' @param size  a numeric vector of length two giving the point and label size, 
#' respectively.
#' @param id.n  an integer giving the number of the most extreme observations 
#' to be identified by a label.  The default is to use the number of identified 
#' outliers, which can be different for the different plots.  See 
#' \dQuote{Details} for more information.
#' @param \dots  for the generic function \code{diagnosticPlot}, additional 
#' arguments to be passed down to methods.  For the \code{"sparseLTSGrid"} 
#' method of \code{diagnosticPlot}, additional arguments to be passed down to 
#' the \code{"sparseLTS"} method.  For the \code{"sparseLTS"} method of 
#' \code{diagnosticPlot}, additional arguments to be passed down to 
#' \code{\link[ggplot2]{geom_point}}.  For the \code{"sparseLTS"} method of 
#' \code{plot}, additional arguments to be passed down to \code{diagnosticPlot}.
#' 
#' @return  
#' If only one plot is requested, an object of class \code{"ggplot"} (see 
#' \code{\link[ggplot2]{ggplot}}), otherwise a list of such objects.
#' 
#' @author Andreas Alfons
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
        ask = (which == "all"), size = c(2, 4), id.n = NULL, ...) {
    # initializations
    fit <- match.arg(fit)
    which <- match.arg(which)
    scale <- switch(fit, reweighted=x$scale, raw=x$raw.scale, 
        both=c(reweighted=x$scale, raw=x$raw.scale))
    if(all(scale <= 0)) {
        stop("plots not available (residual scale equal to 0)")
    }
    data <- fortify(x, fit=fit)  # all information required for plotting
    # call functions for selected plots
    if(which == "all") {
        devAskNewPage(ask)  # ask for new page (if requested)
        # residual Q-Q plot
        tmp <- try(rqqPlot(x, data, fit=fit, size=size, id.n=id.n, ...), 
            silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res <- list(rqq=tmp)
        }
        # residuals vs indices plot
        tmp <- try(residualPlot(x, data, fit=fit, abscissa="index", 
                size=size, id.n=id.n, ...), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res$rindex <- tmp
        }
        # residuals vs fitted plot
        tmp <- try(residualPlot(x, data, fit=fit, abscissa="fitted", 
                size=size, id.n=id.n, ...), silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res$rfit <- tmp
        }
        # regression diagnostic plot
        tmp <- try(rdiagPlot(x, data, fit=fit, size=size, id.n=id.n, ...), 
            silent=TRUE)
        if(inherits(tmp, "try-error")) {
            warn <- gsub("Error in", "In", tmp)
            warning(warn, call.=FALSE)
        } else {
            print(tmp)
            res$rdiag <- tmp
        }
    } else if(which == "rqq") {
        # residual Q-Q plot
        res <- rqqPlot(x, data, fit=fit, size=size, id.n=id.n, ...)
        print(res)
    } else if(which == "rindex") {
        # residuals vs indices plot
        res <- residualPlot(x, data, fit=fit, abscissa="index", 
            size=size, id.n=id.n, ...)
        print(res)
    } else if(which == "rfit") {
        # residuals vs fitted plot
        res <- residualPlot(x, data, fit=fit, abscissa="fitted", 
            size=size, id.n=id.n, ...)
        print(res)
    } else if(which == "rdiag") {
        # regression diagnostic plot
        res <- rdiagPlot(x, data, fit=fit, size=size, id.n=id.n, ...)
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
    x$raw.fitted.values <- fitted(x, s=s, fit="raw")
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

# ----------------------

## residual Q-Q plot

rqqPlot <- function(x, data, ...) UseMethod("rqqPlot")

rqqPlot.sparseLTS <- function(x, data, fit = c("reweighted", "raw", "both"), 
        size = c(2, 4), id.n = NULL, ...) {
    ## initializations
    fit <- match.arg(fit)
    ## extract data frame for reference line
    lineData <- attr(data, "qqLine")
    ## construct data frame for labels
    labelData <- labelify(data, which="qqd", id.n=id.n)
    ## call workhorse function
    p <- ggRqqPlot(data, lineData, labelData, size=size, ...)
    if(fit == "both") {
        f <- as.formula(paste(".", "fit", sep="~"))
        p <- p + facet_grid(f)
    }
    p
}

## workhorse function
ggRqqPlot <- function(data, lineData = NULL, labelData = NULL, 
        size = c(2, 4), ..., mapping, main, xlab, ylab) {
    # initializations
    size <- as.numeric(size)
    size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
    size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
    # define aesthetic mapping for Q-Q plot
    mapping <- aes_string(x="theoretical", y="residual", color="classification")
    # define default axis labels
    if(missing(main)) main <- "Normal Q-Q plot"
    if(missing(xlab)) xlab <- "Quantiles of the standard normal distribution"
    if(missing(ylab)) ylab <- "Standardized sparse LTS residual"
    # create plot
    p <- ggplot(data)
    if(!is.null(lineData)) {
        # add reference line
        lineMapping <- aes_string(intercept="intercept", slope="slope")
        p <- p + geom_abline(lineMapping, lineData, alpha=0.4)
    }
    p <- p + geom_point(mapping, size=size[1], ...) 
    if(!is.null(labelData)) {
        # add labels for observations with largest distances
        labelMapping <- aes_string(x="theoretical", y="residual", label="index")
        p <- p + geom_text(labelMapping, data=labelData, 
            hjust=0, size=size[2], alpha=0.4)
    }
    p <- p + opts(title=main) + labs(x=xlab, y=ylab)
    p
}

## compute theoretical quantiles
qqNorm <- function(y) {
    # TODO: NA handling
    n <- length(y)                # number of observations
    prob <- ppoints(n)            # probabilities
    qnorm(prob)[order(order(y))]  # theoretical quantiles in original order
}

## compute intercept and slope of reference line
qqLine <- function(y) {
    prob <- c(0.25, 0.75)
    ly <- quantile(y, prob, na.rm=TRUE, names=FALSE)
    lx <- qnorm(prob)
    slope <- diff(ly) / diff(lx)
    intercept <- ly[1] - slope * lx[1]
    list(intercept=intercept, slope=slope)
}

# ----------------------

## plot standardized residuals vs indices or fitted values

residualPlot <- function(x, data, ...) UseMethod("residualPlot")

residualPlot.sparseLTS <- function(x, data, 
        fit = c("reweighted", "raw", "both"), 
        abscissa = c("index", "fitted"), 
        size = c(2, 4), id.n = NULL, ...) {
    ## initializations
    fit <- match.arg(fit)
    ## construct data frame for labels
    labelData <- labelify(data, which="residual", id.n=id.n)
    ## call workhorse function
    p <- ggResidualPlot(data, labelData, abscissa=abscissa, size=size, ...)
    if(fit == "both") {
        f <- as.formula(paste(".", "fit", sep="~"))
        p <- p + facet_grid(f)
    }
    p
}

## workhorse function
ggResidualPlot <- function(data, labelData = NULL, 
        abscissa = c("index", "fitted"), size = c(2, 4), 
        ..., mapping, main, xlab, ylab) {
    # initializations
    abscissa <- match.arg(abscissa)
    size <- as.numeric(size)
    size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
    size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
    # define aesthetic mapping for Q-Q plot
    mapping <- aes_string(x=abscissa, y="residual", color="classification")
    # define default axis labels
    if(missing(main)) {
        postfix <- switch(abscissa, index="indices", fitted="fitted values")
        main <- paste("Residuals vs", postfix)
    }
    if(missing(xlab)) {
        xlab <- switch(abscissa, index="Index", fitted="Fitted value")
    }
    if(missing(ylab)) ylab <- "Standardized sparse LTS residual"
    # ensure that horizontal grid line is drawn at 0
    breaks <- union(pretty(data[, "residual"]), 0)
    # create plot
    p <- ggplot(data) + 
        geom_hline(aes(yintercept=-2.5), alpha=0.4) + 
        geom_hline(aes(yintercept=2.5), alpha=0.4) + 
        geom_point(mapping, size=size[1], ...) 
    if(!is.null(labelData)) {
        # add labels for observations with largest distances
        labelMapping <- aes_string(x=abscissa, y="residual", label="index")
        p <- p + geom_text(labelMapping, data=labelData, 
            hjust=0, size=size[2], alpha=0.4)
    }
    p <- p + scale_y_continuous(breaks=breaks) + 
        opts(title=main) + labs(x=xlab, y=ylab)
    p
}

# ----------------------

## plot robust distances

rdiagPlot <- function(x, data, ...) UseMethod("rdiagPlot")

rdiagPlot.sparseLTS <- function(x, data, fit = c("reweighted", "raw", "both"), 
       size = c(2, 4), id.n = NULL, ...) {
    ## initializations
    fit <- match.arg(fit)
    if(fit == "both") 
        onlyNA <- tapply(is.na(data[, "rd"]), data[, "fit", drop=FALSE], all)
    else onlyNA <- all(is.na(data[, "rd"]))
    if(any(onlyNA)) stop("robust distances not available")
    ## extract data frame for vertical reference line
    lineData <- attr(data, "q")
    ## construct data frame for labels
    labelData <- labelify(data, which="xyd", id.n=id.n)
    ## call workhorse function
    p <- ggRdiagPlot(data, lineData, labelData, size=size, ...)
    if(fit == "both") {
        f <- as.formula(paste(".", "fit", sep="~"))
        p <- p + facet_grid(f)
    }
    p
}

## workhorse function
ggRdiagPlot <- function(data, lineData = NULL, labelData = NULL, 
        size = c(2, 4), ..., mapping, main, xlab, ylab) {
    # initializations
    size <- as.numeric(size)
    size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
    size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
    # define aesthetic mapping for Q-Q plot
    mapping <- aes_string(x="rd", y="residual", color="classification")
    # define default axis labels
    if(missing(main)) main <- "Regression diagnostic plot"
    if(missing(xlab)) xlab <- "Robust distance computed by MCD"
    if(missing(ylab)) ylab <- "Standardized sparse LTS residual"
    # create plot
    p <- ggplot(data) + 
        geom_hline(aes(yintercept=-2.5), alpha=0.4) + 
        geom_hline(aes(yintercept=2.5), alpha=0.4)
    if(!is.null(lineData)) {
        # add reference line
        p <- p + geom_vline(aes_string(xintercept="q"), lineData, alpha=0.4)
    }
    p <- p + geom_point(mapping, size=size[1], ...) 
    if(!is.null(labelData)) {
        # add labels for observations with largest distances
        labelMapping <- aes_string(x="rd", y="residual", label="index")
        p <- p + geom_text(labelMapping, data=labelData, 
            hjust=0, size=size[2], alpha=0.4)
    }
    p <- p + opts(title=main) + labs(x=xlab, y=ylab)
    p
}
