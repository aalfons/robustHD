/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_UTILS_H
#define _robustHD_UTILS_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>

using namespace arma;

// functions to export to R
RcppExport SEXP R_findSmallest(SEXP R_x, SEXP R_h);
RcppExport SEXP R_partialOrder(SEXP R_x, SEXP R_h);
//SEXP R_partialSort(SEXP R_x, SEXP R_h);

// functions to be used within C++
vec applyScaleFun(const mat& x, SEXP scaleFun);
uvec findSmallest(const vec& x, const uword& h);
uvec partialOrder(const vec& x, const uword& h);
uvec seqLen(const uword& n);
sword sign(const double& x);

#endif
