# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## take subsets of sparse LTS models over a grid of lambda values

subset.sparseLTSGrid <- function(x, s = NULL, drop = TRUE, ...) {
    # initializations
    if(is.null(s)) return(x)
    call <- match.call()
    call[[1]] <- as.name("subset")
    # extract subset from each component
    x$coefficients <- coef(x, s=s, fit="reweighted", drop=drop)
    x$fitted.values <- fitted(x, s=s, fit="reweighted", drop=drop)
    x$residuals <- residuals(x, s=s, fit="reweighted", drop=drop)
    x$wt <- wt(x, s=s, fit="reweighted", drop=drop)
    x$raw.coefficients <- coef(x, s=s, fit="raw", drop=drop)
    x$raw.fitted.values <- fitted(x, s=s, fit="raw", drop=drop)
    x$raw.residuals <- residuals(x, s=s, fit="raw", drop=drop)
    x$raw.wt <- wt(x, s=s, fit="raw", drop=drop)
    x$best <- x$best[, s, drop=drop]
    x$objective <- x$objective[s]
    x$center <- x$center[s]
    x$scale <- x$scale[s]
    x$lambda <- x$lambda[s]
    x$cnp2 <- x$cnp2[s]
    x$df <- x$df[s]
    x$raw.center <- x$raw.center[s]
    x$raw.scale <- x$raw.scale[s]
    x$raw.cnp2 <- x$raw.cnp2[s]
    if(length(s) == 1 && isTRUE(drop)) {
        # remove names from vectors
        x$objective <- unname(x$objective)
        x$center <- unname(x$center)
        x$scale <- unname(x$scale)
        x$lambda <- unname(x$lambda)
        x$cnp2 <- unname(x$cnp2)
        x$df <- unname(x$df)
        x$raw.center <- unname(x$raw.center)
        x$raw.scale <- unname(x$raw.scale)
        x$raw.cnp2 <- unname(x$raw.cnp2)
        # drop unnecessary components
        x$crit <- x$critValues <- x$sOpt <- x$raw.critValues <- x$raw.sOpt <- NULL
        class(x) <- "sparseLTS"
    } else {
        # select new optimal models
        x$critValues <- x$critValues[s]
        x$raw.critValues <- x$raw.critValues[s]
        if(x$crit == "BIC") {
            x$sOpt <- which.min(x$critValues)
            x$raw.sOpt <- which.min(x$raw.critValues)
        }
    }
    # return subset
    x$call <- call
    x
}
