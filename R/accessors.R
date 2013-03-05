# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## get component 'best'
# generic function
getBest <- function(x, ...) UseMethod("getBest")

# get index of final model according to BIC
getBest.bicSelect <- function(x, ...) x$best

# return NULL by default
getBest.default <- function(x, ...) NULL


## get a component for certain steps of the model sequence
# this is used for accessors that are exported to the namespace 
# (coefficients, fitted values, residuals, ...), so checks for the 
# arguments are necessary

# generic function
getComponent <- function(x, which, ...) UseMethod("getComponent")

# method for class "sparseLTS"
getComponent.sparseLTS <- function(x, which, s = NA, 
                                   fit = c("reweighted", "raw", "both"), 
                                   drop = !is.null(s), ...) {
  ## initializations
  fit <- match.arg(fit)
  if(fit != "reweighted") raw.which <- paste("raw", which, sep=".")
  sMax <- length(x$lambda)
  ## check lambda contains more than one value for the penalty parameter 
  if(sMax > 1) {
    # extract component
    comp <- switch(fit, reweighted=x[[which]], raw=x[[raw.which]], 
                   both={
                     rew <- x[[which]]
                     raw <- x[[raw.which]]
                     if(is.null(dim(rew))) unlist(list(reweighted=rew, raw=raw))
                     else {
                       colnames(rew) <- paste("reweighted", 
                                              colnames(rew), sep=".")
                       colnames(raw) <- paste("raw", colnames(raw), sep=".")
                       cbind(rew, raw)
                     }
                   })
    # check selected steps and extract corresponding parts of the component
    if(!is.null(s)) {
      if(isTRUE(is.na(s))) {
        s <- getSOpt(x, fit=fit)  # optimal step size as default
        if(fit == "both") s[2] <- sMax + s[2]
      } else if(fit == "both" && is.list(s)) {
        # list of steps for each fit
        s <- rep(s, length.out=2)
        s <- lapply(s, checkSteps, sMin=1, sMax=sMax)
        s <- c(s[[1]], sMax+s[[2]])
      } else {
        s <- checkSteps(s, sMin=1, sMax=sMax)
        if(fit == "both") s <- c(s, sMax+s)
      }
      if(is.null(dim(comp))) comp <- comp[s]  # extract selected steps
      else {
        comp <- comp[, s, drop=FALSE]           # extract selected steps
        if(isTRUE(drop)) comp <- dropCol(comp)  # drop dimension if requested
      }
    }
  } else {
    # extract component
    comp <- switch(fit, reweighted=x[[which]], raw=x[[raw.which]],
                   both={
                     rew <- x[[which]]
                     raw <- x[[raw.which]]
                     cfun <- if(length(rew) > 1) cbind else c
                     cfun(reweighted=rew, raw=raw)
                   })
  }
  ## return component
  comp
}


## get residual scale
getScale <- function(x, ...) UseMethod("getScale")

getScale.sparseLTS <- function(x, s = NA, fit = "reweighted", ...) {
  getComponent(x, "scale", s=s, fit=fit, ...)
}


## get optimal step
getSOpt <- function(x, ...) UseMethod("getSOpt")

getSOpt.sparseLTS <- function(x, fit = "reweighted", ...) {
  sOpt <- getBest(x$crit)
  if(fit != "both") sOpt <- sOpt[fit]
  sOpt
}
