/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_FASTRLARS_H
#define _robustHD_FASTRLARS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "utils.h"
#include "corHuber.h"

// functions to export to R
RcppExport SEXP R_fastRlars(SEXP R_x, SEXP R_y, SEXP R_sMax, SEXP R_c,
		SEXP R_prob, SEXP R_tol);

#endif
