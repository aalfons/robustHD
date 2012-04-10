# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## for testing purposes
fastLasso <- function(x, y, lambda, subset = NULL, intercept = TRUE, 
        eps = .Machine$double.eps, use.Gram = TRUE) {
    # initializations
    if(is.null(subset)) {
        useSubset <- FALSE
        subset <- integer()
    } else {
        useSubset <- TRUE
        subset <- subset - 1
    }
    # call C++ function
    callBackend <- getBackend()
    callBackend("R_fastLasso", R_x=x, R_y=y, R_lambda=lambda, 
        R_useSubset=useSubset, R_subset=subset, R_intercept=isTRUE(intercept), 
        R_eps=eps, R_useGram=isTRUE(use.Gram))
}
