/*
 * Author: Andreas Alfons
 *         Erasmus Universiteit Rotterdam
 */

#ifndef _robustHD_INITIALSUBSETS_H
#define _robustHD_INITIALSUBSETS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "fastLasso.h"
#include "utils.h"

// functions to export to R
RcppExport
SEXP R_sparseSubsets(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_h,
                     SEXP R_subsets, SEXP R_normalize, SEXP R_intercept,
                     SEXP R_eps, SEXP R_useGram);

#endif
