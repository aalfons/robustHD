# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Plot a sequence of regression models
#' 
#' Produce a plot of the coefficients or values of the optimality criterion for 
#' a sequence of regression models, such as submodels along a robust least 
#' angle regression sequence, or sparse least trimmed squares regression models 
#' for a grid of values for the penalty parameter.
#' 
#' @method plot seqModel
#' @aliases plot.rlars
#' 
#' @param x  the model fit to be plotted.
#' @param method  a character string specifying the type of plot.  Possible 
#' values are \code{"coefficients"} to plot the coefficients from the submodels 
#' via \code{\link{coefPlot}}, \code{"crit"} to plot the values of the 
#' optimality criterion for the submodels via \code{\link{critPlot}}, or 
#' \code{"diagnostic"} for diagnostic plots via \code{\link{diagnosticPlot}}.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return  
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{coefPlot}}, \code{\link{critPlot}}, 
#' \code{\link{diagnosticPlot}}, \code{\link{rlars}}, \code{\link{sparseLTS}}
#' 
#' @example inst/doc/examples/example-plot.R
#' 
#' @keywords hplot
#' 
#' @export

plot.seqModel <- function(x, method = c("coefficients", "crit", "diagnostic"), 
                          ...) {
  ## initializations
  method <- if(is.null(getSOpt(x))) "diagnostic" else match.arg(method)
  ## call plot function
  if(method == "coefficients") coefPlot(x, ...)
  else if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}


#' @rdname plot.seqModel
#' @method plot sparseLTS
#' @export 

plot.sparseLTS <- function(x, method = c("coefficients", "crit", "diagnostic"), 
                           ...) {
  ## initializations
  crit <- x$crit
  if(is.null(crit) || inherits(crit, "fitSelect")) method <- "diagnostic" 
  else method <- match.arg(method)
  ## call plot function
  if(method == "coefficients") coefPlot(x, ...)
  else if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}

# ----------------------

## supplement the coefficients in a model with other useful information
## returns a data frame suitable for plotting with ggplot2

coefify <- function(model, ...) UseMethod("coefify")

coefify.seqModel <- function(model, zeros = FALSE, labels, ...) {
  # prepare coefficients and labels
  coef <- removeIntercept(t(coef(model, s=NULL)))
  sigmaX <- model$sigmaX
  if(!isTRUE(zeros)) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop=FALSE]
    sigmaX <- sigmaX[keep]
    if(!is.null(labels)) labels <- labels[keep]
  }
  # standardize coefficients
  coef <- sweep(coef, 2, sigmaX, "/", check.margin=FALSE)
  # prepare other information
  m <- ncol(coef)          # number of variables
  steps <- model$s         # step numbers
  nsteps <- length(steps)  # number of steps
  df <- model$df           # degrees of freedom
  vn <- colnames(coef)     # variable names
  # build data frame
  coefData <- data.frame(step=rep.int(steps, m), 
                         df=rep.int(df, m), coefficient=as.numeric(coef), 
                         variable=factor(rep(vn, each=nsteps), levels=vn))
  if(!is.null(labels)) 
    coefData$label <- rep(as.character(labels), each=nsteps)
  coefData
}

coefify.sparseLTS <- function(model, fit = c("reweighted", "raw", "both"), 
                              zeros = FALSE, labels, ...) {
  # initializations
  fit <- match.arg(fit)
  zeros <- isTRUE(zeros)
  coef <- removeIntercept(t(coef(model, s=NULL, fit=fit)))
  df <- getComponent(model, "df", s=NULL, fit=fit)
  # prepare coefficients and labels
  if(!zeros) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop=FALSE]
    if(!is.null(labels)) labels <- labels[keep]
  }
  # check if predictor data is available to compute scale estimates
  if(is.null(x <- model$x)) {
    x <- try(model.matrix(model$terms), silent=TRUE)
    if(inherits(x, "try-error")) 
      stop("scale estimates of predictor variables not available")
  }
  x <- removeIntercept(x)
  if(!zeros) x <- x[, keep, drop=FALSE]
#   # obtain scale estimates for predictors
#   lambda <- model$lambda      # tuning parameters
#   steps <- seq_along(lambda)  # step numbers
#   if(fit %in% c("reweighted", "both")) {
#     cdelta <- model$cnp2
#     wt <- as.matrix(wt(model, s=NULL, fit="raw"))
#     sigmaX <- do.call(rbind, 
#                       lapply(steps, function(s) {
#                         xOk <- x[wt[, s] == 1, , drop=FALSE]
#                         apply(xOk, 2, sd) * cdelta[s]
#                       }))
#   } else sigmaX <- NULL
#   if(fit %in% c("raw", "both")) {
#     cdelta <- model$raw.cnp2
#     best <- as.matrix(model$best)
#     raw.sigmaX <- do.call(rbind, 
#                           lapply(steps, function(s) {
#                             xBest <- x[best[, s], , drop=FALSE]
#                             apply(xBest, 2, sd) * cdelta
#                           }))
#   } else raw.sigmaX <- NULL
#   sigmaX <- rbind(sigmaX, raw.sigmaX)
#   # standardize coeffients
#   coef <- coef / sigmaX
#   # prepare other information
#   m <- ncol(coef)        # number of variables
#   sMax <- length(steps)  # number of steps
#   vn <- colnames(coef)   # variable names
  # obtain scale estimates for predictors
  n <- nrow(x)
  sigmaX <- apply(x, 2, function(x) {
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
  m <- ncol(coef)             # number of variables
  lambda <- model$lambda      # tuning parameters
  steps <- seq_along(lambda)  # step numbers
  sMax <- length(steps)       # number of steps
  vn <- colnames(coef)        # variable names
  # build data frame
  if(fit == "both") {
    fits <- c("reweighted", "raw")
    coefData <- data.frame(
      fit=rep.int(factor(rep(fits, each=sMax), levels=fits), m), 
      lambda=rep.int(lambda, 2*m), step=rep.int(steps, 2*m), 
      df=rep.int(df, m), coefficient=as.numeric(coef), 
      variable=factor(rep(vn, each=2*sMax), levels=vn))
    if(!is.null(labels)) 
      coefData$label <- rep(as.character(labels), each=2*sMax)
  } else {
    coefData <- data.frame(
      lambda=rep.int(lambda, m), step=rep.int(steps, m), 
      df=rep.int(df, m), coefficient=as.numeric(coef), 
      variable=factor(rep(vn, each=sMax), levels=vn))
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
#' @method coefPlot sparseLTS
#' @export

coefPlot.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"), 
                               abscissa = c("step", "df"), zeros = FALSE, 
                               size = c(0.5, 2, 4), labels, offset = 1, ...) {
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
    # split plot into different panels
    p <- p + facet_grid(. ~ fit)
  }
  p
}


## workhorse function
ggCoefPlot <- function(coefData, labelData, abscissa = c("step", "df"), 
                       zeros = FALSE, size = c(0.5, 2, 4), labels, offset = 1, 
                       main = NULL, xlab, ylab, ..., mapping, data) {
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
    theme(legend.position="none") + 
    labs(title=main, x=xlab, y=ylab)
}

# ----------------------
