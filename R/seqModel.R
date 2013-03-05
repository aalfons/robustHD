# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

seqModel <- function(x, y, active, sMin = 0, sMax = NULL, assign = NULL, 
                     robust = TRUE, regFun = .lmrob.fit, useFormula = FALSE, 
                     regArgs = list(), crit = "BIC", cl = NULL) {
  # initializations
  n <- length(y)
  haveAssign <- !is.null(assign)
  if(haveAssign && !is.list(assign)) {
    # column indices for each block in list form
    assign <- split(seq_len(length(assign)), assign)
  }
  if(robust) callRegFun <- getCallFun(regArgs)
#   # add ones to matrix of predictors to account for intercept
#   x <- addIntercept(x, check=TRUE)
  # define steps along the sequence for which to compute submodels
  # the default is to fit models as long as there are twice as many 
  # observations as predictors
  if(haveAssign) {
    if(is.null(sMax)) {
      dfMax <- floor(n/2) + 1
      sMax <- min(length(active), dfMax - 1)
    } else {
      # if 'sMax' is supplied, keep only those steps for which there are 
      # more observations than predictors
      dfMax <- n
    }
  } else if(is.null(sMax)) sMax <- min(length(active), floor(n/2))
  if(sMin > sMax) sMin <- sMax
  s <- sMin:sMax
  # prepare the variable sequence and the degrees of freedom of the models
  if(haveAssign) {
    # compute degrees of freedom of the submodels along sequence
    whichMax <- which.max(s)
    firstActive <- head(active, s[whichMax])
    p <- sapply(assign[firstActive], length)  # number of variables in each block
    df <- unname(cumsum(c(1, p)))[s+1]        # degrees of freedom
    # only fit submodels while the degrees of freedom does not become 
    # larger than the requested maximum
    # FIXME: this could result in an empty sequence
    if(df[whichMax] > dfMax) {
      keep <- which(df <= dfMax)
      s <- s[keep]
      whichMax <- which.max(s)
      firstActive <- head(active, s[whichMax])
      df <- df[keep]
    }
    # groupwise sequenced variables (including intercept)
    sequenced <- c(1, unlist(assign[firstActive], use.names=FALSE) + 1)
  } else {
    # compute degrees of freedom of the submodels along sequence
    df <- s + 1  # account for intercept
    # sequenced variables (including intercept)
    sequenced <- c(1, head(active, max(s)) + 1)
  }
  # define function to fit the submodels along the sequence
  if(robust) {
    if(useFormula) {
      fitFun <- function(k) {
        x <- x[, sequenced[seq_len(k)], drop=FALSE]
        callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
      }
    } else {
      fitFun <- function(k) {
        x <- x[, sequenced[seq_len(k)], drop=FALSE]
        callRegFun(x, y, fun=regFun, args=regArgs)
      }
    }
  } else {
    fitFun <- function(k) lm.fit(x[, sequenced[seq_len(k)], drop=FALSE], y)
  }
  # fit submodels
  # number of variables to use is one less than degrees of freedom
  if(is.null(cl)) models <- lapply(df, fitFun) 
  else models <- parLapply(cl, df, fitFun)
  # construct matrix of coefficents
  coef <- matrix(0, nrow=ncol(x), ncol=length(s), 
                 dimnames=list(colnames(x), s))
  for(k in seq_len(ncol(coef))) {
    coef[sequenced[seq_len(df[k])], k] <- coef(models[[k]])
  }
  # extract fitted values and residuals along sequence
  fitted <- sapply(models, fitted)
  residuals <- sapply(models, residuals)
  colnames(fitted) <- colnames(residuals) <- s
  # construct return object
  out <- list(active=active, s=s, df=df, coefficients=coef, 
              fitted.values=fitted, residuals=residuals, 
              robust=robust)
  class(out) <- "seqModel"
  # add robust scale estimates if applicable
#   if(robust) out$scale <- sapply(models, function(x) getScale(x))
  if(robust) out$scale <- sapply(models, getScale)
  # compute optimality criterion along sequence
  if(crit == "BIC") out$crit <- bicSelect(out)
  # return results
  out
}
