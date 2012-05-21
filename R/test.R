# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## for testing purposes

# test object-based lasso
testLasso <- function(x, y, lambda, initial, intercept = TRUE, 
        eps = .Machine$double.eps, use.Gram = TRUE) {
    callBackend <- getBackend()
    fit <- callBackend("R_testLasso", R_x=x, R_y=y, R_lambda=lambda, 
        R_initial=initial, R_intercept=intercept, R_eps=eps, 
        R_useGram=use.Gram)
    fit$indices <- drop(fit$indices) + 1
    fit$coefficients <- drop(fit$coefficients)
    fit$residuals <- drop(fit$residuals)
    fit
}

# test object-based C-step
testCStep <- function(x, y, lambda, subset, h, intercept = TRUE, 
        tol = .Machine$double.eps^0.5, eps = .Machine$double.eps, 
        use.Gram = TRUE) {
    callBackend <- getBackend()
    fit <- callBackend("R_testCStep", R_x=x, R_y=y, R_lambda=lambda, 
        R_subset=subset, R_intercept=intercept, R_tol=tol, R_eps=eps, 
        R_useGram=use.Gram)
    fit$indices <- drop(fit$indices) + 1
    fit$coefficients <- drop(fit$coefficients)
    fit$residuals <- drop(fit$residuals)
    fit
}

# test keeping best subsets
testKeepBest <- function(subsetMat, crits, nkeep) {
    callBackend <- getBackend()
    out <- callBackend("R_testKeepBest", R_subsetMat=subsetMat, R_crits=crits, 
        R_nkeep=nkeep)
    out$crits <- drop(out$crits)
    out
}
