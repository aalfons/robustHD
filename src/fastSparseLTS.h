/*
 * Author: Andreas Alfons
 *         K.U.Leuven
 */

#ifndef _robustHD_FASTSPARSELTS_H
#define _robustHD_FASTSPARSELTS_H

#include <robustHD.h>

RcppExport SEXP R_findSmallest(SEXP R_x, SEXP R_h);
//RcppExport SEXP R_partialOrder(SEXP R_x, SEXP R_h);
RcppExport SEXP R_initialSubsetsSparse(SEXP R_x, SEXP R_y, SEXP R_subsets,
		SEXP R_h, SEXP R_lambda, SEXP R_intercept, SEXP R_eps, SEXP R_useGram);
RcppExport SEXP R_fastLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_useSubset,
		SEXP R_subset, SEXP R_intercept, SEXP R_eps, SEXP R_useGram);
RcppExport SEXP R_fastCStep(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_residuals,
		SEXP R_h, SEXP R_intercept, SEXP R_eps, SEXP R_useGram);
RcppExport SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subsets,
		SEXP R_intercept, SEXP R_nkeep, SEXP R_ncstep, SEXP R_tol, SEXP R_eps,
		SEXP R_useGram);

#endif
