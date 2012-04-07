/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_CORHUBER_H
#define _robustHD_CORHUBER_H

#define EIGEN_NO_DEBUG

#include <RcppEigen.h>

using namespace Rcpp;

// functions to export to R
RcppExport SEXP R_corHuberUni(SEXP R_x, SEXP R_y, SEXP R_c);
RcppExport SEXP R_corHuberAdj(SEXP R_x, SEXP R_y, SEXP R_c);
RcppExport SEXP R_corHuberBi(SEXP R_x, SEXP R_y, SEXP R_c,
		SEXP R_prob, SEXP R_tol);
RcppExport SEXP R_corMatHuber(SEXP R_x, SEXP R_c, SEXP R_prob, SEXP R_tol);

#endif
