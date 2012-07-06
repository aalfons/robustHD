# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Convert a sparse LTS regression model into a data frame for plotting
#' 
#' Supplement the fitted values and residuals of a sparse least trimmed squared 
#' model with other useful information for diagnostic plots.
#' 
#' @method fortify sparseLTS
#' 
#' @param model  the model fit to be converted.
#' @param data  currently ignored.
#' @param s  an integer vector giving the indices of the models to be 
#' converted.  The default is to use the optimal model for each of the 
#' requested fits.
#' @param fit  a character string specifying which fit to convert.  Possible 
#' values are \code{"reweighted"} (the default) to convert the reweighted fit, 
#' \code{"raw"} to convert the raw fit, or \code{"both"} to convert both fits.
#' @param \dots  currently ignored.
#' 
#' @return  A data frame containing the columns listed below, as well as 
#' additional information stored in the attributes \code{"qqLine"} (intercepts 
#' and slopes of the respective reference lines to be displayed in residual Q-Q 
#' plots), \code{"q"} (quantiles of the Mahalanobis distribution used as cutoff 
#' points for detecting leverage points), and \code{"facets"} (default faceting 
#' formulas for the diagnostic plots).
#' @returnItem step  the indices of the models (only returned if more than one 
#' model is requested).
#' @returnItem fit  the model fits (only returned if both the reweighted 
#' and raw fit are requested).
#' @returnItem index  the indices of the observations.
#' @returnItem fitted  the fitted values.
#' @returnItem residual  the standardized residuals.
#' @returnItem theoretical  the corresponding theoretical quantiles from the 
#' standard normal distribution.
#' @returnItem qqd  the absolute distances from a reference line through the 
#' first and third sample and theoretical quartiles.
#' @returnItem rd  the robust Mahalanobis distances computed via the MCD (see 
#' \code{\link[robustbase]{covMcd}}).
#' @returnItem xyd  the pairwise maxima of the absolute values of the 
#' standardized residuals and the robust Mahalanobis distances, divided by the 
#' respective other outlier detection cutoff point.
#' @returnItem weight  the weights indicating regression outliers.
#' @returnItem leverage  logicals indicating leverage points (i.e., outliers in 
#' the predictor space).
#' @returnItem classification  a factor with levels \code{"outlier"} 
#' (regression outliers) and \code{"good"} (data points following the model).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[ggplot2]{fortify}}, \code{\link{diagnosticPlot}}, 
#' \code{\link{sparseLTS}}, \code{\link{sparseLTSGrid}}
#' 
#' @example inst/doc/examples/example-fortify.R
#' 
#' @keywords utilities
#' 
#' @import ggplot2
#' @export

fortify.sparseLTS <- function(model, data, 
        fit = c("reweighted", "raw", "both"), ...) {
    ## initializations
    fit <- match.arg(fit)
    scale <- switch(fit, reweighted=model$scale, raw=model$raw.scale, 
        both=c(reweighted=model$scale, raw=model$raw.scale))
    if(any(scale <= 0)) stop("residual scale equal to 0")
    ## construct data frame with all information for plotting
    if(fit == "both") {
        fits <- c("reweighted", "raw")
        ## recursive call for each fit
        reweighted <- fortify(model, fit="reweighted")
        raw <- fortify(model, fit="raw")
        ## combine data for Q-Q reference line
        qql <- data.frame(fit=factor(fits, levels=fits), 
            rbind(attr(reweighted, "qqLine"), attr(raw, "qqLine")), 
            row.names=NULL)
        ## combine data for cutoff chi-squared quantile
        q <- data.frame(fit=factor(fits, levels=fits), 
            rbind(attr(reweighted, "q"), attr(raw, "q")), 
            row.names=NULL)
        ## combine results
        n <- c(nrow(reweighted), nrow(raw))
        data <- data.frame(fit=rep.int(factor(fits, levels=fits), n), 
            rbind(reweighted, raw), row.names=NULL)
        attr(data, "facets") <- . ~ fit
        attr(data, "qqLine") <- qql
        attr(data, "q") <- q
    } else {
        ## extract fitted values
        fitted <- fitted(model, fit=fit)
        ## extract standardized residuals
        residuals <- residuals(model, fit=fit, standardized=TRUE)
        n <- length(residuals)  # number of observations
        ## extract outlier weights
        wt <- wt(model, fit=fit)
        ## compute theoretical quantiles and distances from Q-Q reference line
        theoretical <- qqNorm(residuals)
        qql <- qqLine(residuals)  # Q-Q reference line
        qqd <- abs(residuals - qql$intercept - qql$slope * theoretical)
        ## compute MCD distances using significant variables
        # extract predictor matrix
        terms <- delete.response(model$terms)  # extract terms for model matrix
        ok <- TRUE
        if(is.null(x <- model$x)) {
            x <- try(model.matrix(terms), silent=TRUE)
            if(inherits(x, "try-error")) {
                ok <- FALSE
                warning("model data not available")
            }
        }
        if(ok) {
            if(model$intercept) x <- removeIntercept(x)
            # extract coefficients
            coefficients <- coef(model, fit=fit)
            if(model$intercept) coefficients <- removeIntercept(coefficients)
            significant <- which(coefficients != 0)
            p <- length(significant)
            if(p == 0) {
                ok <- FALSE
                warning("all coefficients equal to 0")
            }
        }
        if(ok) {
            # adjust alpha since MCD computes subset size depending on n and p
            h <- model$quan
            n2 <- (n+p+1) %/% 2
            alpha <- pmin((h - 2*n2 + n) / (2 * (n - n2)), 1)
            # check fraction for subset size
            if(alpha < 0.5) {
                alpha <- 0.5
                warning(sprintf("cannot compute MCD with h = %d; using h = %d", 
                        model$quan, h.alpha.n(alpha, n, p)))
            }
            # compute distances
            rd <- try({
                    xs <- x[, significant, drop=FALSE]
                    mcd <- covMcd(xs, alpha=alpha)
                    if(fit == "reweighted") {
                        center <- mcd$center
                        cov <- mcd$cov
                    } else {
                        center <- mcd$raw.center
                        cov <- mcd$raw.cov
                    }
                    sqrt(mahalanobis(xs, center, cov))
                })
            if(inherits(rd, "try-eror")) {
                ok <- FALSE
                warning("robust distances cannot be computed")
            }
        }
        if(!ok) rd <- rep.int(NA, n)
        # take maximum of the distances in the x- and y-space, divided by the 
        # respective other cutoff point
        q <- sqrt(qchisq(0.975, p))
        xyd <- pmax.int(abs(rd/2.5), abs(residuals/q))
        ## construct indicator variables for leverage points
        leverage <- rd > q
        ## classify data points
        class <- ifelse(wt == 0, "outlier", "good")
        class <- factor(class, levels=c("outlier", "good"))
        ## construct data frame
        data <- data.frame(index=seq_len(n), fitted=fitted, residual=residuals, 
            theoretical=theoretical, qqd=qqd, rd=rd, xyd=xyd, weight=wt, 
            leverage=leverage, classification=class)
        attr(data, "qqLine") <- as.data.frame(qql)
        attr(data, "q") <- data.frame(q=max(q, 2.5))
    }
    ## return data frame
    data
}


#' @rdname fortify.sparseLTS
#' @method fortify sparseLTSGrid
#' @export 

fortify.sparseLTSGrid <- function(model, data, s, 
        fit = c("reweighted", "raw", "both"), ...) {
    ## initializations
    fit <- match.arg(fit)
    steps <- getSteps(model)
    lambda <- model$lambda
    # check the steps and scale estimate
    bothOpt <- FALSE
    if(missing(s)) {
        if(fit == "both") {
            s <- unique(c(model$sOpt, model$raw.sOpt))
            bothOpt <- length(s) == 2
        } else s <- switch(fit, reweighted=model$sOpt, raw=model$raw.sOpt)
    } else if(is.null(s)) s <- steps 
    else {
        sMax <- length(steps)
        s <- checkSteps(s, sMin=1, sMax=sMax)
    }
    # check scale estimate
    if(bothOpt) {
        scale <- c(model$scale[s[1]], model$raw.scale[s[2]])
    } else {
        scale <- switch(fit, reweighted=model$scale[s], raw=model$raw.scale[s], 
            both=c(model$scale[s], model$raw.scale[s]))
    }
    if(any(scale <= 0)) stop("residual scale equal to 0")
    ## extract data for the requested steps
    if(bothOpt) {
        # extract the data from the respecitve optimal lambda
        fits <- c("reweighted", "raw")
        ## recursive call for each fit
        reweighted <- fortify(subset(model, s=s[1]), fit="reweighted")
        raw <- fortify(subset(model, s=s[2]), fit="raw")
        ## combine data for Q-Q reference line
        qql <- data.frame(fit=factor(fits, levels=fits), 
            rbind(attr(reweighted, "qqLine"), attr(raw, "qqLine")), 
            row.names=NULL)
        ## combine data for cutoff chi-squared quantile
        q <- data.frame(fit=factor(fits, levels=fits), 
            rbind(attr(reweighted, "q"), attr(raw, "q")), 
            row.names=NULL)
        ## combine results
        n <- c(nrow(reweighted), nrow(raw))
        data <- data.frame(fit=rep.int(factor(fits, levels=fits), n), 
            rbind(reweighted, raw), row.names=NULL)
        attr(data, "facets") <- . ~ fit
        attr(data, "qqLine") <- qql
        attr(data, "q") <- q
    } else if(length(s) == 1) {
        # extract the data from the selected step
        data <- fortify(subset(model, s=s), fit=fit)
    } else {
        # extract the data from each requested step
        data <- lapply(s, function(s) fortify(subset(model, s=s), fit=fit))
        qql <- lapply(data, attr, which="qqLine")
        q <- lapply(data, attr, which="q")
        # combine data from the steps
        data <- cbind(step=rep.int(s, sapply(data, nrow)), do.call(rbind, data))
        qql <- cbind(step=rep.int(s, sapply(qql, nrow)), do.call(rbind, qql))
        q <- cbind(step=rep.int(s, sapply(q, nrow)), do.call(rbind, q))
        attr(data, "facets") <- if(fit == "both") step ~ fit else ~step
        attr(data, "qqLine") <- qql
        attr(data, "q") <- q
    }
    ## return data frame
    data
}
