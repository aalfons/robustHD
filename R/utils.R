# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## add default column names to matrix
addColnames <- function(x) {
  # 'x' needs to be a matrix
  if(is.null(colnames(x))) colnames(x) <- paste("x", seq_len(ncol(x)), sep="")
  x
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
  if(!check || is.na(match("(Intercept)", colnames(x)))) {
    cbind("(Intercept)"=rep.int(1, nrow(x)), x)
  } else x
}

## copy names from a vector or matrix to another vector or matrix
copyNames <- function(from, to, which = "col", target = "row") {
  # read names from source
  if(is.null(dim(from))) nam <- names(from) 
  else if(which == "row") nam <- rownames(from)
  else if(which == "col") nam <- colnames(from)
  # write names to target
  if(is.null(dim(to))) names(to) <- nam
  else if(target == "row") rownames(to) <- nam
  else if(target == "col") colnames(to) <- nam
  # return object
  to
}

## find indices of h smallest observations
findSmallest <- function(x, h) {
  # call C++ function
  callBackend <- getBackend()
  callBackend("R_findSmallest", R_x=as.numeric(x), R_h=as.integer(h))
}

## compute coefficients of hyperplane through data points
hyperplane <- function(x) {
  p <- ncol(x)
  y <- -x[, p]  # right hand side
  x <- x[, -p, drop=FALSE]
  tx <- t(x)
  theta <- solve(tx %*% x) %*% tx %*% y
  c(theta, 1)
}

## obtain the degrees of freedom of a model (number of nonzero parameters)
modelDf <- function(beta, tol = .Machine$double.eps^0.5) {
  length(which(abs(beta) > tol))
}
