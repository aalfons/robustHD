/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _robustHD_FASTSPARSELTS_H
#define _robustHD_FASTSPARSELTS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>	// OpenMP
#endif
#include "fastLasso.h"
#include "utils.h"

// functions to export to R
RcppExport SEXP R_testLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
  	SEXP R_normalize, SEXP R_intercept, SEXP R_eps, SEXP R_useGram);
RcppExport SEXP R_testCStep(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subset,
  	SEXP R_normalize, SEXP R_intercept, SEXP R_tol, SEXP R_eps, 
    SEXP R_useGram);
RcppExport SEXP R_testKeepBest(SEXP R_subsetMat, SEXP R_crits, SEXP R_nkeep);
RcppExport SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
  	SEXP R_normalize, SEXP R_intercept, SEXP R_ncstep, SEXP R_nkeep, 
    SEXP R_tol, SEXP R_eps, SEXP R_useGram, SEXP R_ncores);

#endif
