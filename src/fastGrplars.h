/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_FASTGRPLARS_H
#define _robustHD_FASTGRPLARS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>	// OpenMP
#endif
#include "utils.h"
#include "corHuber.h"

// functions to export to R
RcppExport SEXP R_fastGrplars(SEXP R_x, SEXP R_y, SEXP R_sMax,
		SEXP R_assign, SEXP R_ncores);

#endif
