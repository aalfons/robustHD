# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## for testing purposes
fastCStep <- function(x, y, lambda, residuals, h, intercept = TRUE, 
        eps = .Machine$double.eps, useGram = TRUE) {
    fit <- .Call("R_fastCStep", R_x=x, R_y=y, R_lambda=lambda, 
        R_residuals=residuals, R_h=h, R_intercept=intercept, R_eps=eps, 
        R_useGram=useGram, PACKAGE="robustHD")
    fit$subset <- fit$subset + 1
    fit
}
