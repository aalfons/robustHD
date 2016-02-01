/*
 * Authors: Andreas Alfons
 *          Erasmus Universiteit Rotterdam
 *
 *          Viktoria Oellerer
 *          KU Leuven
 */

#ifndef _robustHD_FASTROBUSTLASSO_H
#define _robustHD_FASTROBUSTLASSO_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>	// OpenMP
#endif
#include "fastLasso.h"
#include "Mscale.h"
#include "utils.h"

// functions to export to R
RcppExport
SEXP R_rawLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_normalize,
                SEXP R_intercept, SEXP R_eps, SEXP R_useGram);
RcppExport
SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
                     SEXP R_normalize, SEXP R_intercept, SEXP R_ncstep,
                     SEXP R_nkeep, SEXP R_tol, SEXP R_eps, SEXP R_useGram,
                     SEXP R_ncores);
RcppExport
SEXP R_fastSparseS(SEXP R_x, SEXP R_y, SEXP R_findLambda, SEXP R_fixedLambda,
                   SEXP R_initial, SEXP R_nistep, SEXP R_nkeep, SEXP R_k,
                   SEXP R_b, SEXP R_nfpi, SEXP R_nfpiMax, SEXP R_tol,
                   SEXP R_eps, SEXP R_useGram, SEXP R_ncores);
RcppExport
SEXP R_fastSparseMM(SEXP R_x, SEXP R_y, SEXP R_findLambda, SEXP R_fixedLambda,
                    SEXP R_fixedLambdaS, SEXP R_k, SEXP R_rMax, SEXP R_initial,
                    SEXP R_nistep, SEXP R_nkeep, SEXP R_kS, SEXP R_b,
                    SEXP R_nfpi, SEXP R_nfpiMax, SEXP R_tol, SEXP R_eps,
                    SEXP R_useGram, SEXP R_ncores);

#endif
