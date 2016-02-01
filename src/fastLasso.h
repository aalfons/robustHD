/*
 * Author: Andreas Alfons
 *         Erasmus Universiteit Rotterdam
 */

#ifndef _robustHD_FASTLASSO_H
#define _robustHD_FASTLASSO_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "utils.h"

using namespace arma;

// functions to export to R
RcppExport
SEXP R_fastLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_useSubset,
                 SEXP R_subset, SEXP R_normalize, SEXP R_intercept,
                 SEXP R_eps, SEXP R_useGram);
RcppExport
SEXP R_testLasso(SEXP R_x, SEXP R_y, SEXP R_findLambda, SEXP R_lambda,
                 SEXP R_useWeights, SEXP R_weights, SEXP R_tol,
                 SEXP R_eps, SEXP R_useGram);

// functions to be used within C++
void fastLasso(const mat& x, const vec& y, const double& lambda,
               const bool& useSubset, const uvec& subset, const bool& normalize,
               const bool& useIntercept, const double& eps, const bool& useGram,
               const bool& useCrit,
               // intercept, coefficients, residuals and objective function are
               // returned through the following parameters
               double& intercept, vec& beta, vec& residuals, double& crit);
void fastLasso(const mat& x, const vec& y, const bool& findLambda,
               double& lambda, const bool& useWeights, const vec& weights,
               const double& tol, const double& eps, const bool& useGram,
               // intercept, coefficients and residuals are returned through
               // the following parameters
               double& intercept, vec& beta, vec& residuals);

#endif
