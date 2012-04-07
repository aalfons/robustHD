/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include "initialSubsets.h"

using namespace Rcpp;
using namespace Eigen;

// function used internally by initialSubsets(), which computes lasso fits for
// subsets containing a small number of observations (typically only 3) and
// returns the indices of the respective h observations with the smallest
// absolute residuals
MatrixXi initialSubsetsSparse(const MatrixXd& x, const VectorXd& y,
		const MatrixXi& subsets, const int& h, const double& lambda,
		const bool& useIntercept, const double& eps,
		const bool& useGram) {
	const int nsamp = subsets.cols();
	MatrixXi indices(h, nsamp);
	for(int k = 0; k < nsamp; k++) {
		// compute lasso fit
		double intercept;
		VectorXd coefficients = fastLasso(x, y, lambda, true, subsets.col(k),
				useIntercept, eps, useGram, intercept);
		// compute residuals
		VectorXd residuals;
		residuals.noalias() = y - x * coefficients;
		if(useIntercept) {
			for(int i = 0; i < residuals.size(); i++) {
				residuals(i) -= intercept;
			}
		}
		// find h observations with smallest absolute residuals
		indices.col(k) = findSmallest(residuals.cwiseAbs(), h);
	}
	return indices;
}

// R interface to initialSubsetsSparse()
SEXP R_initialSubsetsSparse(SEXP R_x, SEXP R_y, SEXP R_subsets, SEXP R_h,
		SEXP R_lambda, SEXP R_intercept, SEXP R_eps, SEXP R_useGram) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	IntegerMatrix Rcpp_subsets(R_subsets);	// subset to use for computation
	const int nsamp = Rcpp_subsets.ncol();
	Map<MatrixXi> subsets(Rcpp_subsets.begin(), Rcpp_subsets.nrow(), nsamp);
	int h = as<int>(R_h);
	double lambda = as<double>(R_lambda);
	bool useIntercept = as<bool>(R_intercept);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function and return results
	MatrixXi indices = initialSubsetsSparse(x, y, subsets, h, lambda,
			useIntercept, eps, useGram);
	return wrap(indices);
}
