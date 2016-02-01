/*
 * Author: Andreas Alfons
 *         Erasmus Universiteit Rotterdam
 */

#ifndef _robustHD_FASTLARS_H
#define _robustHD_FASTLARS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>	// OpenMP
#endif
#include "utils.h"
#include "corHuber.h"

// functions to export to R
RcppExport SEXP R_fastLars(SEXP R_x, SEXP R_y, SEXP R_sMax, SEXP R_robust, SEXP R_c,
		SEXP R_prob, SEXP R_tol, SEXP scaleFun, SEXP R_ncores);

#endif
