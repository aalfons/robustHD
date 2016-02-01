# ---------------------------------------
# Authors: Viktoria Oellerer
#          KU Leuven
#
#          Andreas Alfons
#          Erasmus Universiteit Rotterdam
# ---------------------------------------

# biweight rho-function for M-estimator
Mrho <- function(x, k = 1.547645) {
  # call C++ function
  out <- .Call("R_Mrho", R_x=x, R_k=k, PACKAGE="robustHD")
  # return values
  drop(out)
}

# biweight weight function for M-estimator
Mwgt <- function(x, k = 1.547645) {
  # call C++ function
  w <- .Call("R_Mwgt", R_x=x, R_k=k, PACKAGE="robustHD")
  # return values
  drop(w)
}

# M-estimator of scale
Mscale <- function(x, k = 1.547645, b = NULL, nfpi = 500,
                   tol = .Machine$double.eps^0.5, warn = TRUE) {
  # call C++ function
  if(is.null(b)) b <- Tbsb(k, 1)/k * Mrho(k, k)
  scale <- .Call("R_Mscale", R_x=x, R_k=k, R_b=b, R_nfpi=nfpi,
                 R_tol=tol, R_warn=warn, PACKAGE="robustHD")
  # return values
  drop(scale)
}
