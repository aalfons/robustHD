# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

# lasso on subset of observations
fastLasso <- function(x, y, lambda, subset = NULL, normalize = TRUE,
                      intercept = TRUE, eps = .Machine$double.eps,
                      use.Gram = TRUE, raw = FALSE, drop = TRUE) {
  # initializations
  normalize <- isTRUE(normalize)
  intercept <- isTRUE(intercept)
  use.Gram <- isTRUE(use.Gram)
  raw <- isTRUE(raw)
  drop <- isTRUE(drop)
  # compute lasso
  if(raw) {
    # call C++ function
    fit <- callBackend("R_rawLasso", R_x=x, R_y=y, R_lambda=lambda,
                       R_normalize=normalize, R_intercept=intercept,
                       R_eps=eps, R_useGram=use.Gram)
  } else {
    # check subset
    if(is.null(subset)) {
      useSubset <- FALSE
      subset <- integer()
    } else useSubset <- TRUE
    # call C++ function
    fit <- callBackend("R_subsetLasso", R_x=x, R_y=y, R_lambda=lambda,
                       R_useSubset=useSubset, R_subset=subset,
                       R_normalize=normalize, R_intercept=intercept,
                       R_eps=eps, R_useGram=use.Gram)
  }
  # return lasso fit
  if(drop) fit <- lapply(fit, drop)  # drop the dimension of the components
  fit
}

# weighted lasso (data are assumed to be standardized)
weightedLasso <- function(x, y, weights = NULL, lambdaMin = 0,
                          sMax = NA, tol = .Machine$double.eps^2,
                          eps = .Machine$double.eps, use.Gram = TRUE,
                          drop = TRUE) {
  # initializations
  use.Gram <- isTRUE(use.Gram)
  drop <- isTRUE(drop)
  # check weights
  if(is.null(weights)) {
    n <- length(y)
    weights <- rep.int(1, n)
  } else useWeights <- TRUE
  # check maximum number of steps
  if(is.na(sMax)) sMax <- ncol(x)
  # call C++ function
  fit <- callBackend("R_weightedLasso",
                     R_x=x, R_y=y, R_useWeights=TRUE, R_weights=weights,
                     R_lambdaMin=lambdaMin, R_sMax=sMax, R_tol=tol, R_eps=eps,
                     R_useGram=use.Gram)
  # return lasso fit
  if(drop) fit <- lapply(fit, drop)  # drop the dimension of the components
  fit
}
