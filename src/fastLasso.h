/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _robustHD_FASTLASSO_H
#define _robustHD_FASTLASSO_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "utils.h"

using namespace arma;

// functions to export to R
RcppExport SEXP R_fastLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_useSubset,
  	SEXP R_subset, SEXP R_normalize, SEXP R_intercept, SEXP R_eps, 
    SEXP R_useGram);

// functions to be used within C++
vec fastLasso(const mat& x, const vec& y, const double& lambda,
  	const bool& useSubset, const uvec& subset, const bool& normalize, 
    const bool& useIntercept, const double& eps, const bool& useGram, 
    double& intercept);

#endif
