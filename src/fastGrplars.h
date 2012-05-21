/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_FASTGRPLARS_H
#define _robustHD_FASTGRPLARS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "utils.h"

// functions to export to R
RcppExport SEXP R_fastGrplars(SEXP R_x, SEXP R_y, SEXP R_sMax, SEXP R_assign);

#endif
