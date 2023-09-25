# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Plot a sequence of regression models
#'
#' Produce a plot of the coefficients, the values of the optimality criterion,
#' or diagnostic plots for a sequence of regression models, such as submodels
#' along a robust or groupwise least angle regression sequence, or sparse least
#' trimmed squares regression models for a grid of values for the penalty
#' parameter.
#'
#' @method plot seqModel
#' @aliases plot.rlars plot.grplars plot.tslarsP
#'
#' @param x  the model fit to be plotted.
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param method  a character string specifying the type of plot.  Possible
#' values are \code{"coefficients"} to plot the coefficients from the submodels
#' via \code{\link{coefPlot}} (only for the \code{"seqModel"} and
#' \code{"sparseLTS"} methods), \code{"crit"} to plot the values of the
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
#' \code{\link{diagnosticPlot}}, \code{\link{rlars}}, \code{\link{grplars}},
#' \code{\link{rgrplars}}, \code{\link{tslarsP}}, \code{\link{rtslarsP}},
#' \code{\link{tslars}}, \code{\link{rtslars}}, \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-plot.R
#'
#' @keywords hplot
#'
#' @import ggplot2
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
#' @method plot perrySeqModel
#' @export

plot.perrySeqModel <- function(x, method = c("crit", "diagnostic"), ...) {
  ## initializations
  method <- match.arg(method)
  ## call plot function
  if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}


#' @rdname plot.seqModel
#' @method plot tslars
#' @export

plot.tslars <- function(x, p, method = c("coefficients", "crit", "diagnostic"),
                        ...) {
  ## initializations
  method <- match.arg(method)
  ## call plot function
  if(method == "coefficients") coefPlot(x, p=p, ...)
  else if(method == "crit") critPlot(x, p=p, ...)
  else diagnosticPlot(x, p=p, ...)
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


#' @rdname plot.seqModel
#' @method plot perrySparseLTS
#' @export

plot.perrySparseLTS <- function(x, method = c("crit", "diagnostic"), ...) {
  ## initializations
  method <- match.arg(method)
  ## call plot function
  if(method == "crit") critPlot(x, ...)
  else diagnosticPlot(x, ...)
}

# ----------------------


#' Coefficient plot of a sequence of regression models
#'
#' Produce a plot of the coefficients from a sequence of regression models,
#' such as submodels along a robust or groupwise least angle regression
#' sequence, or sparse least trimmed squares regression models for a grid of
#' values for the penalty parameter.
#'
#' @aliases coefPlot.rlars coefPlot.grplars coefPlot.tslarsP
#'
#' @param object  the model fit to be plotted.
#' @param abscissa  a character string specifying what to plot on the
#' \eqn{x}-axis.  For objects inheriting from class \code{"seqModel"}, possible
#' values are \code{"step"} for the step number (the default), or \code{"df"}
#' for the degrees of freedom.  For code{"sparseLTS"} objects, possible values
#' are code{"lambda"} for the value of the penalty parameter (the default), or
#' \code{"step"} for the step number.
#' @param zeros  a logical indicating whether predictors that never enter the
#' model and thus have zero coefficients should be included in the plot
#' (\code{TRUE}) or omitted (\code{FALSE}, the default).  This is useful if the
#' number of predictors is much larger than the number of observations, in
#' which case many coefficients are never nonzero.
#' @param size  a numeric vector of length three giving the line width, the
#' point size and the label size, respectively.
#' @param labels  an optional character vector containing labels for the
#' predictors.  Plotting labels can be suppressed by setting this to
#' \code{NA}.
#' @param offset   an integer giving the offset of the labels from the
#' corresponding coefficient values from the last step (i.e., the number of
#' blank characters to be prepended to the label).
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying for which estimator to produce the
#' plot.  Possible values are \code{"reweighted"} (the default) for the
#' reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for both
#' estimators.
#' @param facets  a faceting formula to override the default behavior.  If
#' supplied, \code{\link[ggplot2]{facet_wrap}} or
#' \code{\link[ggplot2]{facet_grid}} is called depending on whether the formula
#' is one-sided or two-sided.
#' @param \dots  additional arguments to be passed down, eventually to
#' \code{\link[ggplot2]{geom_line}} and \code{\link[ggplot2]{geom_point}}.
#'
#' @return
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}},
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
#' \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-coefPlot.R
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @importFrom utils head tail
#' @export

coefPlot <- function(object, ...) UseMethod("coefPlot")


#' @rdname coefPlot
#' @method coefPlot seqModel
#' @export

coefPlot.seqModel <- function(object, zeros = FALSE, labels = NULL, ...) {
  # extract all information required for plotting
  setup <- setupCoefPlot(object, zeros = zeros, labels = labels)
  # call method for object with all information required for plotting
  coefPlot(setup, ...)
}


#' @rdname coefPlot
#' @method coefPlot tslars
#' @export

coefPlot.tslars <- function(object, p, zeros = FALSE, labels = NULL, ...) {
  # extract all information required for plotting
  if (missing(p)) setup <- setupCoefPlot(object, zeros = zeros, labels = labels)
  else setup <- setupCoefPlot(object, p = p, zeros = zeros, labels = labels)
  # call method for object with all information required for plotting
  coefPlot(setup, ...)
}


#' @rdname coefPlot
#' @method coefPlot sparseLTS
#' @export

coefPlot.sparseLTS <- function(object, fit = c("reweighted", "raw", "both"),
                               zeros = FALSE, labels = NULL, ...) {
  # extract all information required for plotting
  setup <- setupCoefPlot(object, fit = fit, zeros = zeros, labels = labels)
  # call method for object with all information required for plotting
  coefPlot(setup, ...)
}


#' @rdname coefPlot
#' @method coefPlot setupCoefPlot
#' @export

coefPlot.setupCoefPlot <- function(object, abscissa = NULL,
                                   size = c(0.5, 2, 4), offset = 1,
                                   facets = object$facets, ...) {
  # initializations
  if (is.null(abscissa)) abscissa <- object$abscissa[1]
  else abscissa <- match.arg(abscissa, choices = object$abscissa)
  size <- as.numeric(size)
  size <- c(size, rep.int(NA, max(0, 3-length(size))))[1:3]  # ensure length 3
  size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
  # define default axis labels
  xlab <- switch(abscissa, step = "Step", lambda = "lambda",
                 df = "Degrees of freedom")
  ylab <- "Standardized coefficients"
  # define aesthetic mapping for plotting coefficients
  coefMapping <- aes_string(x = abscissa, y = "coefficient",
                            color = "variable")
  # define aesthetic mapping for plotting x-axis grid and labels
  if (object$includeLabels) {
    offset <- paste(rep.int(" ", offset), collapse = "")  # whitespace
    labelData <- object$labels
    labelData$label <- paste(offset, labelData$label, sep = "")
    labelMapping <- aes_string(x = abscissa, y = "coefficient",
                               label = "label")
  }
  # create plot
  p <- ggplot() +
    local_geom_line(coefMapping, data = object$coefficients,
                    size = size[1], ...) +
    local_geom_point(coefMapping, data = object$coefficients,
                     size = size[2], ...)
  if (object$includeLabels) {
    p <- p +
      geom_text(labelMapping, data = labelData, hjust = 0,
                size = size[3], alpha = 0.4)
  }
  p <- p +
    labs(x = xlab, y = ylab)
  if (abscissa == "lambda") p <- p + scale_x_reverse()
  if (!is.null(facets)) {
    # split plot into different panels
    if (length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  # return plot
  p
}


# local geoms to override defaults

local_geom_line <- function(..., show.legend = FALSE) {
  geom_line(..., show.legend = show.legend)
}

local_geom_point <- function(..., show.legend = FALSE) {
  geom_point(..., show.legend = show.legend)
}


# ----------------------


#' Optimality criterion plot of a sequence of regression models
#'
#' Produce a plot of the values of the optimality criterion for a sequence of
#' regression models, such as submodels along a robust or groupwise least angle
#' regression sequence, or sparse least trimmed squares regression models for
#' a grid of values for the penalty parameter.
#'
#' @aliases critPlot.rlars critPlot.grplars critPlot.tslarsP
#'
#' @param object  the model fit to be plotted, , or an object containing
#' all necessary information for plotting (as generated by
#' \code{\link{setupCritPlot}}).
#' @param which  a character string specifying the type of plot.  Possible
#' values are \code{"line"} (the default) to plot the (average) results for
#' each model as a connected line, \code{"dot"} to create a dot plot,
#' \code{"box"} to create a box plot, or \code{"density"} to create a smooth
#' density plot.  Note that the last two plots are only available in case of
#' prediction error estimation via repeated resampling.
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying for which estimator to produce the
#' plot.  Possible values are \code{"reweighted"} (the default) for the
#' reweighted fits, \code{"raw"} for the raw fits, or \code{"both"} for both
#' estimators.
#' @param \dots  additional arguments to be passed down, eventually to
#' \code{\link[ggplot2]{geom_line}}, \code{\link[ggplot2]{geom_pointrange}},
#' \code{\link[ggplot2]{geom_boxplot}}, or \code{\link[ggplot2]{geom_density}}.
#'
#' @return
#' An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
#'
#' @note Function \code{\link[perry]{perryPlot}} is used to create the plot,
#' even if the optimality criterion does not correspond to resampling-based p
#' rediction error estimation.  While this can be seen as as a misuse of its
#' functionality, it ensures that all optimality criteria are displayed in the
#' same way.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[perry]{perryPlot}},
#' \code{\link{rlars}}, \code{\link{grplars}}, \code{\link{rgrplars}},
#' \code{\link{tslarsP}}, \code{\link{rtslarsP}}, \code{\link{tslars}},
#' \code{\link{rtslars}}, \code{\link{sparseLTS}}
#'
#' @example inst/doc/examples/example-critPlot.R
#'
#' @keywords hplot
#'
#' @export

critPlot <- function(object, ...) UseMethod("critPlot")


#' @rdname critPlot
#' @method critPlot seqModel
#' @export

critPlot.seqModel <- function(object, which = c("line", "dot"), ...) {
  # extract all information required for plotting
  setup <- setupCritPlot(object, which = which)
  # call method for object with all information required for plotting
  critPlot(setup, ...)
}


#' @rdname critPlot
#' @method critPlot tslars
#' @export

critPlot.tslars <- function(object, p, which = c("line", "dot"), ...) {
  # extract all information required for plotting
  if (missing(p)) setup <- setupCritPlot(object, which = which)
  else setup <- setupCritPlot(object, p = p, which = which)
  # call method for object with all information required for plotting
  critPlot(setup, ...)
}


#' @rdname critPlot
#' @method critPlot sparseLTS
#' @export

critPlot.sparseLTS <- function(object, which = c("line", "dot"),
                               fit = c("reweighted", "raw", "both"),
                               ...) {
  # extract all information required for plotting
  setup <- setupCritPlot(object, which = which, fit = fit)
  # call method for object with all information required for plotting
  critPlot(setup, ...)
}


#' @rdname critPlot
#' @method critPlot perrySeqModel
#' @export

critPlot.perrySeqModel <- function(object,
                                   which = c("line", "dot", "box", "density"),
                                   ...) {
  # extract all information required for plotting
  setup <- setupCritPlot(object, which = which)
  # call method for object with all information required for plotting
  critPlot(setup, ...)
}


#' @rdname critPlot
#' @method critPlot perrySparseLTS
#' @export

critPlot.perrySparseLTS <- function(object,
                                    which = c("line", "dot", "box", "density"),
                                    fit = c("reweighted", "raw", "both"),
                                    ...) {
  # extract all information required for plotting
  setup <- setupCritPlot(object, which = which, fit = fit)
  # call method for object with all information required for plotting
  critPlot(setup, ...)
}


#' @rdname critPlot
#' @method critPlot setupCritPlot
#' @export

critPlot.setupCritPlot <- function(object, ...) {
  # define x-axis label
  tuning <- object$tuning
  useTuning <- !is.null(tuning) && ncol(tuning) == 1
  xlab <- if (useTuning) names(tuning) else "Step"
  # make sure y-axis label is correct
  if (inherits(object, "setupBICPlot")) ylab <- "BIC"
  else if (inherits(object, "setupPerryPlot")) ylab <- "Prediction error"
  else ylab <- "Optimality criterion"  # shouldn't happen
  # initialize plot (this is a bit of a hack)
  class(object) <- "setupPerryPlot"
  p <- perryPlot(object, ...) + labs(x = xlab, y = ylab)
  # reverse x-axis for penalty parameter of sparse LTS
  if (xlab == "lambda") p <- p + scale_x_reverse()
  # return plot
  p
}


# ----------------------


## construct data frame for labels based on some order
getLabelData <- function(data, which, id.n = NULL) {
  # initializations
  if(isTRUE(id.n < 1)) return(NULL)
  by <- intersect(c("step", "fit"), names(data))
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


#' Diagnostic plots for a sequence of regression models
#'
#' Produce diagnostic plots for a sequence of regression models, such as
#' submodels along a robust least angle regression sequence, or sparse least
#' trimmed squares regression models for a grid of values for the penalty
#' parameter.  Four plots are currently implemented.
#'
#' In the normal Q-Q plot of the standardized residuals, a reference line is
#' drawn through the first and third quartile.  The \code{id.n} observations
#' with the largest distances from that line are identified by a label (the
#' observation number).  The default for \code{id.n} is the number of
#' regression outliers, i.e., the number of observations whose residuals are
#' too large (cf. \code{\link[=weights.sparseLTS]{weights}}).
#'
#' In the plots of the standardized residuals versus their index or the fitted
#' values, horizontal reference lines are drawn at 0 and +/-2.5.  The
#' \code{id.n} observations with the largest absolute values of the
#' standardized residuals are identified by a label (the observation
#' number).  The default for \code{id.n} is the number of regression outliers,
#' i.e., the number of observations whose absolute residuals are too large (cf.
#' \code{\link[=weights.sparseLTS]{weights}}).
#'
#' For the regression diagnostic plot, the robust Mahalanobis distances of the
#' predictor variables are computed via the minimum covariance determinant
#' (MCD) estimator based on only those predictors with non-zero coefficients
#' (see \code{\link[robustbase]{covMcd}}).  Horizontal reference lines are
#' drawn at +/-2.5 and a vertical reference line is drawn at the upper 97.5\%
#' quantile of the \eqn{\chi^{2}}{chi-squared} distribution with \eqn{p}
#' degrees of freedom, where \eqn{p} denotes the number of predictors with
#' non-zero coefficients.  The \code{id.n} observations with the largest
#' absolute values of the standardized residuals and/or largest robust
#' Mahalanobis distances are identified by a label (the observation number).
#' The default for \code{id.n} is the number of all outliers: regression
#' outliers (i.e., observations whose absolute residuals are too large, cf.
#' \code{\link[=weights.sparseLTS]{weights}}) and leverage points (i.e.,
#' observations with robust Mahalanobis distance larger than the 97.5\%
#' quantile of the \eqn{\chi^{2}}{chi-squared} distribution with \eqn{p}
#' degrees of freedom).
#'
#' Note that the argument \code{alpha} for controlling the subset size
#' behaves differently for \code{\link{sparseLTS}} than for
#' \code{\link[robustbase]{covMcd}}.  For \code{\link{sparseLTS}}, the subset
#' size \eqn{h} is determined by the fraction \code{alpha} of the number of
#' observations \eqn{n}.  For \code{\link[robustbase]{covMcd}}, on the other
#' hand, the subset size also depends on the number of variables \eqn{p} (see
#' \code{\link[robustbase]{h.alpha.n}}).  However, the \code{"sparseLTS"} and
#' \code{"perrySparseLTS"} methods attempt to compute the MCD using the same
#' subset size that is used to compute the sparse least trimmed squares
#' regressions.  This may not be possible if the number of selected variables
#' is large compared to the number of observations. In such cases,
#' \code{\link{setupDiagnosticPlot}} returns \code{NA}s for the robust
#' Mahalanobis distances, and the regression diagnostic plot fails.
#'
#' @aliases diagnosticPlot.rlars diagnosticPlot.grplars diagnosticPlot.tslarsP
#'
#' @param object  the model fit for which to produce diagnostic plots, or an
#' object containing all necessary information for plotting (as generated
#' by \code{\link{setupDiagnosticPlot}}).
#' @param s  for the \code{"seqModel"} method, an integer vector giving
#' the steps of the submodels  for which to produce diagnostic plots (the
#' default is to use the optimal submodel).  For the \code{"sparseLTS"} method,
#' an integer vector giving the indices of the models for which to produce
#' diagnostic plots (the default is to use the optimal model for each of the
#' requested fits).
#' @param covArgs  a list of arguments to be passed to
#' \code{\link[robustbase]{covMcd}} for the regression diagnostic plot (see
#' @param p  an integer giving the lag length for which to produce the plot
#' (the default is to use the optimal lag length).
#' @param fit  a character string specifying for which fit to produce
#' diagnostic plots.  Possible values are \code{"reweighted"} (the default) for
#' diagnostic plots for the reweighted fit, \code{"raw"} for diagnostic plots
#' for the raw fit, or \code{"both"} for diagnostic plots for both fits.
#' \dQuote{Details}).
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
#' @param facets  a faceting formula to override the default behavior.  If
#' supplied, \code{\link[ggplot2]{facet_wrap}} or
#' \code{\link[ggplot2]{facet_grid}} is called depending on whether the formula
#' is one-sided or two-sided.
#' @param size  a numeric vector of length two giving the point and label size,
#' respectively.
#' @param id.n  an integer giving the number of the most extreme observations
#' to be identified by a label.  The default is to use the number of identified
#' outliers, which can be different for the different plots.  See
#' \dQuote{Details} for more information.
#' @param \dots  additional arguments to be passed down, eventually to
#' \code{\link[ggplot2]{geom_point}}.
#'
#' @return
#' If only one plot is requested, an object of class \code{"ggplot"} (see
#' \code{\link[ggplot2]{ggplot}}), otherwise a list of such objects.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link{rlars}},
#' \code{\link{grplars}}, \code{\link{rgrplars}}, \code{\link{tslarsP}},
#' \code{\link{rtslarsP}}, \code{\link{tslars}}, \code{\link{rtslars}},
#' \code{\link{sparseLTS}}, \code{\link[robustbase:ltsPlot]{plot.lts}}
#'
#' @example inst/doc/examples/example-diagnosticPlot.R
#'
#' @keywords hplot
#'
#' @import ggplot2
#' @import robustbase
#' @importFrom grDevices devAskNewPage
#' @export

diagnosticPlot <- function(object, ...) UseMethod("diagnosticPlot")


#' @rdname diagnosticPlot
#' @method diagnosticPlot seqModel
#' @export

diagnosticPlot.seqModel <- function(object, s = NA, covArgs = list(), ...) {
  # extract all information required for plotting
  setup <- setupDiagnosticPlot(object, s = s, covArgs = covArgs)
  # call method for object with all information required for plotting
  diagnosticPlot(setup, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot perrySeqModel
#' @export

diagnosticPlot.perrySeqModel <- function(object, covArgs = list(), ...) {
  # extract all information required for plotting
  setup <- setupDiagnosticPlot(object, covArgs = covArgs)
  # call method for object with all information required for plotting
  diagnosticPlot(setup, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot tslars
#' @export

diagnosticPlot.tslars <- function(object, p, s = NA, covArgs = list(), ...) {
  # extract all information required for plotting
  if (missing(p)) setup <- setupDiagnosticPlot(object, s = s, covArgs = covArgs)
  else setup <- setupDiagnosticPlot(object, p = p, s = s, covArgs = covArgs)
  # call method for object with all information required for plotting
  diagnosticPlot(setup, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot sparseLTS
#' @export

diagnosticPlot.sparseLTS <- function(object, s = NA,
                                     fit = c("reweighted", "raw", "both"),
                                     covArgs = list(), ...) {
  # extract all information required for plotting
  setup <- setupDiagnosticPlot(object, s = s, fit = fit, covArgs = covArgs)
  # call method for object with all information required for plotting
  diagnosticPlot(setup, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot perrySparseLTS
#' @export

diagnosticPlot.perrySparseLTS <- function(object,
                                          fit = c("reweighted", "raw", "both"),
                                          covArgs = list(), ...) {
  # extract all information required for plotting
  setup <- setupDiagnosticPlot(object, fit = fit, covArgs = covArgs)
  # call method for object with all information required for plotting
  diagnosticPlot(setup, ...)
}


#' @rdname diagnosticPlot
#' @method diagnosticPlot setupDiagnosticPlot
#' @export

diagnosticPlot.setupDiagnosticPlot <- function(object,
                                               which = c("all", "rqq",
                                                         "rindex", "rfit",
                                                         "rdiag"),
                                               ask = (which == "all"),
                                               facets = object$facets,
                                               size = c(2, 4), id.n = NULL,
                                               ...) {
  # initializations
  which <- match.arg(which)
  size <- as.numeric(size)
  size <- c(size, rep.int(NA, max(0, 2-length(size))))[1:2]  # ensure length 2
  size <- ifelse(is.na(size), eval(formals()$size), size)    # fill NA's
  # call functions for selected plots
  if (which == "all") {
    oldAsk <- devAskNewPage(ask)  # ask for new page (if requested)
    on.exit(devAskNewPage(oldAsk))
    # residual Q-Q plot
    p <- try(rqqPlot(object, facets = facets, size = size, id.n = id.n, ...),
             silent = TRUE)
    if (inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call. = FALSE)
      res <- list()
    } else {
      print(p)
      res <- list(rqq = p)
    }
    # residuals vs indices plot
    p <- try(residualPlot(object, abscissa = "index", facets = facets,
                          size = size, id.n = id.n, ...),
             silent = TRUE)
    if (inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call. = FALSE)
    } else {
      print(p)
      res$rindex <- p
    }
    # residuals vs fitted plot
    p <- try(residualPlot(object, abscissa = "fitted", facets = facets,
                          size = size, id.n = id.n, ...),
             silent = TRUE)
    if (inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call. = FALSE)
    } else {
      print(p)
      res$rfit <- p
    }
    # regression diagnostic plot
    p <- try(rdiagPlot(object, facets = facets, size = size, id.n = id.n, ...),
             silent = TRUE)
    if (inherits(p, "try-error")) {
      warn <- gsub("Error in", "In", p)
      warning(warn, call. = FALSE)
    } else {
      print(p)
      res$rdiag <- p
    }
    invisible(res)
  } else if (which == "rqq") {
    # residual Q-Q plot
    rqqPlot(object, facets = facets, size = size, id.n = id.n, ...)
  } else if (which == "rindex") {
    # residuals vs indices plot
    residualPlot(object, abscissa = "index", facets = facets,
                 size = size, id.n = id.n, ...)
  } else if (which == "rfit") {
    # residuals vs fitted plot
    residualPlot(object, abscissa = "fitted", facets = facets,
                 size = size, id.n = id.n, ...)
  } else if (which == "rdiag") {
    # regression diagnostic plot
    rdiagPlot(object, facets = facets, size = size, id.n = id.n, ...)
  }
}

# ----------------------

rqqPlot <- function(object, facets = object$facets, size = c(2, 4),
                    id.n = NULL, ..., mapping) {
  # define aesthetic mapping for Q-Q plot
  mapping <- aes_string(x = "theoretical", y = "residual",
                        color = "Diagnostics")
  # extract data frame for reference line
  lineData <- object$qqLine
  # construct data frame for labels
  labelData <- getLabelData(object$data, which = "qqd", id.n = id.n)
  # define default title and axis labels
  main <- "Normal Q-Q plot"
  xlab <- "Quantiles of the standard normal distribution"
  ylab <- "Standardized residual"
  # create plot
  p <- ggplot(object$data)
  if (!is.null(lineData)) {
    # add reference line
    lineMapping <- aes_string(intercept = "intercept", slope = "slope")
    p <- p + geom_abline(lineMapping, lineData, alpha = 0.4)
  }
  p <- p + geom_point(mapping, size = size[1], ...)
  if(!is.null(labelData)) {
    # add labels for observations with largest distances
    labelMapping <- aes_string(x = "theoretical", y = "residual",
                               label = "index")
    p <- p + geom_text(labelMapping, data = labelData, hjust = 0,
                       size = size[2], alpha = 0.4)
  }
  p <- p + labs(title = main, x = xlab, y = ylab)
  if (!is.null(facets)) {
    # split plot into different panels
    if (length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
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

residualPlot <- function(object, abscissa = c("index", "fitted"),
                         facets = object$facets, size = c(2, 4),
                         id.n = NULL, ..., mapping) {
  ## initializations
  abscissa <- match.arg(abscissa)
  # define aesthetic mapping for residual plot
  mapping <- aes_string(x = abscissa, y = "residual", color = "Diagnostics")
  ## construct data frame for labels
  labelData <- getLabelData(object$data, which = "residual", id.n = id.n)
  # define default title and axis labels
  postfix <- switch(abscissa, index = "indices", fitted = "fitted values")
  main <- paste("Residuals vs", postfix)
  xlab <- switch(abscissa, index = "Index", fitted = "Fitted value")
  ylab <- "Standardized residual"
  # ensure that horizontal grid line is drawn at 0
  breaks <- union(pretty(object$data[, "residual"]), 0)
  # create plot
  p <- ggplot(object$data) +
    geom_hline(aes(yintercept = -2.5), alpha = 0.4) +
    geom_hline(aes(yintercept = 2.5), alpha = 0.4) +
    geom_point(mapping, size = size[1], ...)
  if (!is.null(labelData)) {
    # add labels for observations with largest distances
    labelMapping <- aes_string(x = abscissa, y = "residual", label = "index")
    p <- p + geom_text(labelMapping, data = labelData, hjust = 0,
                       size = size[2], alpha = 0.4)
  }
  p <- p + scale_y_continuous(breaks = breaks) +
    labs(title = main, x = xlab, y = ylab)
  if (!is.null(facets)) {
    # split plot into different panels
    if (length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  p
}

# ----------------------

## plot robust distances (regression diagnostic plot)

rdiagPlot <- function(object, facets = object$facets, size = c(2, 4),
                      id.n = NULL, ..., mapping) {
  ## initializations
  # extract data frame with main information
  data <- object$data
  # extract data frame for vertical reference line
  lineData <- object$q
  # check if robust distances are available
  msg <- "robust distances not available"
  by <- intersect(c("step", "fit"), names(data))
  if (length(by) > 0) {
    indices <- split(seq_len(nrow(data)), data[, by, drop = FALSE])
    onlyNA <- sapply(indices, function(i) all(is.na(data[i, "rd"])))
    if (all(onlyNA)) stop(msg)
    if (any(onlyNA)) {
      indices <- do.call(c, unname(indices[onlyNA]))
      data <- data[-indices, , drop = FALSE]
      lineData <- lineData[!onlyNA, , drop = FALSE]
      warning(msg, " for some submodels")
    }
  } else {
    onlyNA <- all(is.na(data[, "rd"]))
    if (onlyNA) stop(msg)
  }
  # define aesthetic mapping for regression diagnostic plot
  mapping <- aes_string(x = "rd", y = "residual", color = "Diagnostics")
  ## construct data frame for labels
  labelData <- getLabelData(data, which = "xyd", id.n = id.n)
  # define default title and axis labels
  main <- "Regression diagnostic plot"
  xlab <- "Robust distance computed by MCD"
  ylab <- "Standardized residual"
  # create plot
  p <- ggplot(data) +
    geom_hline(aes(yintercept = -2.5), alpha = 0.4) +
    geom_hline(aes(yintercept = 2.5), alpha = 0.4)
  if(!is.null(lineData)) {
    # add reference line
    p <- p + geom_vline(aes_string(xintercept = "q"), lineData, alpha = 0.4)
  }
  p <- p + geom_point(mapping, size = size[1], ...)
  if (!is.null(labelData)) {
    # add labels for observations with largest distances
    labelMapping <- aes_string(x = "rd", y = "residual", label = "index")
    p <- p + geom_text(labelMapping, data = labelData, hjust = 0,
                       size = size[2], alpha = 0.4)
  }
  p <- p + labs(title = main, x = xlab, y =  ylab)
  if (!is.null(facets)) {
    # split plot into different panels
    if (length(facets) == 2) p <- p + facet_wrap(facets)
    else p <- p + facet_grid(facets)
  }
  p
}
