# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

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
    fit <- callBackend("R_fastLasso", R_x=x, R_y=y, R_lambda=lambda,
                       R_useSubset=useSubset, R_subset=subset,
                       R_normalize=normalize, R_intercept=intercept,
                       R_eps=eps, R_useGram=use.Gram)
  }
  # return lasso fit
  if(drop) fit <- lapply(fit, drop)  # drop the dimension of the components
  fit
}
