/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_FASTSPARSELTS_H
#define _robustHD_FASTSPARSELTS_H

#define EIGEN_NO_DEBUG

#include <robustHD.h>
#include "fastLasso.h"
#include "utils.h"

// functions to export to R
RcppExport SEXP R_testLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subset,
		SEXP R_intercept, SEXP R_eps, SEXP R_useGram);
RcppExport SEXP R_testCStep(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_lassoSubset,
		SEXP R_intercept, SEXP R_tol, SEXP R_eps, SEXP R_useGram);
RcppExport SEXP R_testKeepBest(SEXP R_subsetMat, SEXP R_crits, SEXP R_nkeep);
RcppExport SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subsets,
		SEXP R_intercept, SEXP R_ncstep, SEXP R_nkeep, SEXP R_tol, SEXP R_eps,
		SEXP R_useGram);

#endif
