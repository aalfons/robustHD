# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

testLasso <- function(x, y, lambda = NULL, weights = NULL,
                      tol = .Machine$double.eps^2, eps = .Machine$double.eps,
                      use.Gram = TRUE, drop = TRUE) {
  # initializations
  use.Gram <- isTRUE(use.Gram)
  drop <- isTRUE(drop)
  # check penalty parameter
  if(is.null(lambda)) {
    findLambda <- TRUE
    subset <- integer()
  } else findLambda <- FALSE
  # check weights
  if(is.null(weights)) {
    useWeights <- FALSE
    weights <- numeric()
  } else useWeights <- TRUE
  # call C++ function
  fit <- callBackend("R_testLasso",
                     R_x=x, R_y=y, R_findLambda=findLambda, R_lambda=lambda,
                     R_useWeights=useWeights, R_weights=weights, R_tol=tol,
                     R_eps=eps, R_useGram=use.Gram)
  # return lasso fit
  if(drop) fit <- lapply(fit, drop)  # drop the dimension of the components
  fit
}
