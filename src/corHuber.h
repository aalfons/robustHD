/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_CORHUBER_H
#define _robustHD_CORHUBER_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace arma;

// functions to export to R
RcppExport SEXP R_corHuberUni(SEXP R_x, SEXP R_y, SEXP R_c);
RcppExport SEXP R_corHuberAdj(SEXP R_x, SEXP R_y, SEXP R_c);
RcppExport SEXP R_corHuberBi(SEXP R_x, SEXP R_y, SEXP R_c,
		SEXP R_prob, SEXP R_tol);
RcppExport SEXP R_corMatHuber(SEXP R_x, SEXP R_c, SEXP R_prob, SEXP R_tol);

// functions to be used within C++
double corPearson(const vec& x, const vec& y);
double corHuberBi(const vec& x, const vec& y, const double& c,
		const double& prob, const double& tol);

#endif
