/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _robustHD_FASTLASSO_H
#define _robustHD_FASTLASSO_H

#define EIGEN_NO_DEBUG

#include <robustHD.h>

using namespace Eigen;

// functions to export to R
RcppExport SEXP R_fastLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_useSubset,
		SEXP R_subset, SEXP R_intercept, SEXP R_eps, SEXP R_useGram);

// functions to be used within C++
VectorXd fastLasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useSubset, const VectorXi& subset, const bool& useIntercept,
		const double& eps, const bool& useGram, double& intercept);

#endif
