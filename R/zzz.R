# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## register built-in RcppArmadillo back end when package is loaded
#' @import Rcpp 
#' @import RcppArmadillo
#' @useDynLib robustHD
.onLoad <- function(libname, pkgname) registerBackend()
