# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------


#' @export
setupCoefPlot <- function(object, ...) UseMethod("setupCoefPlot")

#' @export
setupCoefPlot.seqModel <- function(object, zeros = FALSE, labels = NULL, ...) {
  # initializations
  zeros <- isTRUE(zeros)
  # extract and reformat relevant coefficiencts
  coef <- removeIntercept(t(coef(object, s = NULL)))
  m <- ncol(coef)          # number of variables
  # check labels
  includeLabels <- !isTRUE(is.na(labels))
  if (includeLabels) {
    if (is.null(labels)) labels <- defaultLabels(object)
    else if (!(is.numeric(labels) || is.character(labels)) || length(labels) != m) {
      stop(sprintf("labels should be a character vector of length %d", m))
    }
  }
  # prepare coefficients and labels
  sigmaX <- object$sigmaX
  if (!zeros) {
    keep <- apply(coef != 0, 2, any)
    coef <- coef[, keep, drop = FALSE]
    sigmaX <- sigmaX[keep]
    if (includeLabels) labels <- labels[keep]
  }
  # standardize coefficients
  coef <- sweep(coef, 2, sigmaX, "/", check.margin = FALSE)
  # prepare other information
  m <- ncol(coef)          # number of variables
  vn <- colnames(coef)     # variable names
  steps <- object$s        # step numbers
  nsteps <- length(steps)  # number of steps
  df <- object$df          # degrees of freedom
  # build data frame for coefficients
  coefData <- data.frame(step = rep.int(steps, m),
                         df = rep.int(df, m), coefficient = as.numeric(coef),
                         variable = factor(rep(vn, each = nsteps), levels = vn))
  # construct return object
  out <- list(coefficients = coefData, steps = steps, df = df,
              includeLabels = includeLabels)
  # add data frame for labels if requested
  if(includeLabels) {
    out$labels <- data.frame(step = rep.int(steps[nsteps], m),
                             df = rep.int(df[nsteps], m),
                             coefficient = coef[nsteps, ],
                             label = labels)
  }
  # add class and return object
  class(out) <- "setupCoefPlot"
  out
}
