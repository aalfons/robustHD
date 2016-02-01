/*
 * Authors: Viktoria Oellerer
 *          KU Leuven
 *
 *          Andreas Alfons
 *          Erasmus Universiteit Rotterdam
 */

#ifndef _robustHD_MSCALE_H
#define _robustHD_MSCALE_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace arma;

// functions to be used within C++
//double Mrho(const double& x, const double& k);
vec Mrho(const vec& x, const double& k);
//double Mwgt(const double& x, const double& k);
vec Mwgt(const vec& x, const double& k);
double Mscale(const vec& x, const double& k, const double& b, const uword& nfpi,
              const double& tol, const bool& warn);

// functions to export to R
RcppExport
SEXP R_Mrho(SEXP R_x, SEXP R_k);
RcppExport
SEXP R_Mwgt(SEXP R_x, SEXP R_k);
RcppExport
SEXP R_Mscale(SEXP R_x, SEXP R_k, SEXP R_b, SEXP R_nfpi, SEXP R_tol,
              SEXP R_warn);

#endif
