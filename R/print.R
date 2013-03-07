# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------


#' @S3method print bicSelect
print.bicSelect <- function(x, best = TRUE, ...) {
  # print BIC values
  cat("\nBIC:\n")
  print(x$values, ...)
  # print optimal model if requested
  if(isTRUE(best)) {
    best <- x$best
    text <- "Index of best model:"
    if(length(best) == 1) {
      best <- as.matrix(best)
      dimnames(best) <- list(text, "")
    } else cat("\n", text, "\n", sep="")
    print(best, ...)
  }
  # return object invisibly
  invisible(x)
}

#' @S3method print seqModel
print.seqModel <- function(x, zeros = FALSE, ...) {
  # print function call
  if(!is.null(call <- x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }
  # print predictor sequence
  active <- t(x$active)
  steps <- seq_len(ncol(active))
  dimnames(active) <- list("Var", steps)
  cat("\nSequence of moves:\n")
  print(active, ...)
  # print coefficients of optimal submodel
  cat("\nCoefficients of optimal submodel:\n")
  print(coef(x, zeros=zeros), ...)
  # print index of optimal submodel
  cat(sprintf("\nOptimal step: %d\n", getSOpt(x)))
  # return object invisibly
  invisible(x)
}

#' @S3method print sparseLTS
print.sparseLTS <- function(x, fit = c("reweighted", "raw", "both"), 
                            zeros = FALSE, ...) {
  # initializations
  fit <- match.arg(fit)
  lambda <- x$lambda
  sOpt <- getSOpt(x, fit=fit)
  # print function call
  if(!is.null(call <- x$call)) {
    cat("\nCall:\n")
    dput(x$call)
  }
  # print coefficients of sparse LTS model
  coefficients <- coef(x, fit=fit, zeros=zeros)
  if(length(lambda) == 1) cat("\nCoefficients:\n")
  else {
    if(fit == "both") colnames(coefficients) <- c("reweighted", "raw")
    cat("\nCoefficients of optimal submodel:\n")
  }
  print(coefficients, ...)
  # print penalty parameter and robust scale estimate
  scale <- getScale(x, fit=fit)
  if(length(lambda) == 1) {
    text <- c("Penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
      lambda <- as.matrix(lambda)
      dimnames(lambda) <- list(text[1], "")
      print(lambda, ...)
      dimnames(scale) <- list(text[2], colnames(coefficients))
      print(scale, ...)
    } else {
      info <- matrix(c(lambda, scale), 2, 1)
      dimnames(info) <- list(text, "")
      print(info, ...)
    }
  } else {
    lambda <- lambda[sOpt]
    info <- rbind(lambda, scale)
    text <- c("Optimal penalty parameter:", "Residual scale estimate:")
    if(fit == "both") {
      dimnames(info) <- list(text, colnames(coefficients))
      cat("\n")
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
