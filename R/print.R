# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' @S3method print sparseLTS
print.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"), 
                            zeros = FALSE, ...) {
  # initializations
  fit <- match.arg(fit)
  sOpt <- getSOpt(x, fit=fit)
  # print function call
  if(!is.null(call <- x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }
  # print coefficients of sparse LTS model
  coefficients <- coef(x, fit=fit, zeros=zeros)
  if(is.null(sOpt)) cat("\nCoefficients:\n")
  else {
    if(fit == "both") colnames(coefficients) <- c("reweighted", "raw")
    cat("\nCoefficients of optimal submodel:\n")
  }
  print(coefficients, ...)
  # print penalty parameter and robust scale estimate
  lambda <- x$lambda
  scale <- getScale(x, fit=fit)
  if(length(lambda) == 1) {
    text <- c("Penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
      lambda <- as.matrix(lambda)
      dimnames(lambda) <- list(text[1], "")
      print(lambda, ...)
      cat("\n", text[2], "\n", sep="")
      print(scale, ...)
    } else {
      info <- matrix(c(lambda, scale), 2, 1)
      dimnames(info) <- list(text, "")
      print(info, ...)
    }
  } else if(is.null(sOpt)) {
    text <- c("Penalty parameters:", "Residual scale estimates:")
    cat("\n")
    if(fit == "both") {
      cat(text[1], "\n", sep="")
      print(lambda, ...)
      cat("\n", text[2], "\n", sep="")
      scale <- cbind(reweighted=x$scale, raw=x$raw.scale)
      print(scale, ...)
    } else {
      info <- rbind(lambda, scale)
      rownames(info) <- text
      print(info, ...)
    }
  } else {
    lambda <- lambda[sOpt]
    info <- rbind(lambda, scale)
    text <- c("Optimal penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
      cat("\n")
      dimnames(info) <- list(text, colnames(coefficients))
    } else dimnames(info) <- list(text, "")
    print(info, ...)
  }
  # return object invisibly
  invisible(x)
}

#' @S3method print perrySparseLTS
#' @import perry
print.perrySparseLTS <- function(x, ...) {
  # print prediction error results
  perry:::print.perryTuning(x, best=FALSE, ...)
  # print optimal value for penalty parameter
  optimalLambda <- x$tuning[x$best, "lambda"]
  if(length(optimalLambda) > 1) {
    names(optimalLambda) <- names(x$best)
    cat("\nOptimal lambda:\n")
  } else {
    optimalLambda <- matrix(optimalLambda, 
                            dimnames=list("Optimal lambda:", ""))
  }
  print(optimalLambda, ...)
  # return object invisibly
  invisible(x)
}
