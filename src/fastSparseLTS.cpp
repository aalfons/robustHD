/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include "fastSparseLTS.h"

using namespace Rcpp;
using namespace Eigen;
using namespace std;


// *****************
// class definitions
// *****************

class Subset {
public:
	// information on the subset
	VectorXd subset;
	double intercept;
	VectorXd coefficients;
	VectorXd residuals;
	double crit;
	bool continueCSteps;

	// constructors
	Subset();
	Subset(VectorXi&);
	// compute lasso solution and residuals
	void lasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
			const bool& useIntercept, const double& eps, const bool& useGram);
	// compute objective function
	void objective(const double& lambda);
	// perform C-Step
	void cStep(const MatrixXd& x, const VectorXd& y, const double& lambda,
			const bool& useIntercept, const double& eps, const bool& useGram);
};

// constructors
inline Subset::Subset() {
	continueCSteps = true;
}
inline Subset::Subset(VectorXi& subset) {
	// TODO
}

// compute lasso solution and residuals
void Subset::lasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useIntercept, const double& eps, const bool& useGram) {
	// TODO
}

// compute objective function
void objective(const double& lambda) {
	// TODO
}

// perform C-Step
void Subset::cStep(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useIntercept, const double& eps, const bool& useGram) {
	// TODO
}


// ****************************
// sparse least trimmed squares
// ****************************

// sparse LTS objective function (L1 penalized trimmed sum of squared residuals)
double objective(const VectorXd& residuals, const VectorXi& subset,
		const double& lambda, const VectorXd& coefficients) {
	// compute sum of squared residuals for subset
	const int h = subset.size();
	double crit = 0;
	for(int i = 0; i < h; i++) {
		crit += pow(residuals(subset(i)), 2);
	}
	// add L1 penalty on coefficients
	crit += h * lambda * coefficients.lpNorm<1>();
	return crit;
}

// wrapper function used internally by fastSparseLTS(), which returns residuals
// and value of objective function through corresponding parameters
VectorXd fastLasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const VectorXi& subset, const bool& useIntercept, const double& eps,
		const bool& useGram, double& intercept, VectorXd& residuals,
		double& crit) {
	// compute coefficients
	VectorXd coefficients = fastLasso(x, y, lambda, true, subset,
			useIntercept, eps, useGram, intercept);
	// compute residuals
	residuals.noalias() = y - x * coefficients;
	if(useIntercept) {
		for(int i = 0; i < residuals.size(); i++) {
			residuals(i) -= intercept;
		}
	}
	// compute value of objective function
	crit = objective(residuals, subset, lambda, coefficients);
	return coefficients;
}

// C-step for sparse LTS for k-th subset
// Eigen library is used for linear algebra
// ATTENTION: residuals, subset, intercept and value of objective function are
//            returned through corresponding parameters
// x .............. predictor matrix
// y .............. response
// lambda ......... penalty parameter
// residuals ...... current residuals (updated through this parameter)
// h .............. subset size
// subset ......... current subset (updated through this parameter)
// useIntercept ... logical indicating whether intercept should be included
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
// intercept ...... intercept is returned through this parameter
// crit ........... value of objective function is returned through this parameter
VectorXd fastCStep(const MatrixXd& x, const VectorXd& y, const double& lambda,
		VectorXd& residuals, const int& h, VectorXi& subset,
		const bool& useIntercept, const double& eps, const bool& useGram,
		double& intercept, double& crit) {
	// update subset
	subset = findSmallest(residuals.cwiseAbs(), h);
	// compute lasso solution for new subset
	return fastLasso(x, y, lambda, subset, useIntercept,
			eps, useGram, intercept, residuals, crit);
}

// R interface to fastCStep() (for testing purposes)
SEXP R_fastCStep(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_residuals,
		SEXP R_h, SEXP R_intercept, SEXP R_eps, SEXP R_useGram) {
    // data initializations
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	double lambda = as<double>(R_lambda);
	NumericVector Rcpp_residuals(R_residuals);
	VectorXd residuals(n);
	for(int i = 0; i < n; i++) {
		residuals(i) = Rcpp_residuals[i];	// do not reuse memory
	}
	int h = as<int>(R_h);
	bool useIntercept = as<bool>(R_intercept);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function and return results as list
	VectorXi subset;
	double intercept;
	double crit;
	VectorXd coefficients = fastCStep(x, y, lambda, residuals, h, subset,
			useIntercept, eps, useGram, intercept, crit);
	NumericVector R_coefficients = wrap(coefficients);
	if(useIntercept) {
		R_coefficients.push_front(intercept);	// prepend intercept
	}
	return List::create(
			Named("coefficients") = R_coefficients,
			Named("subset") = subset,
			Named("residuals") = residuals,
			Named("crit") = crit
			);
}





// obtain the indices of the nkeep best subsets among a matrix of subsets
VectorXi bestSubsets(const VectorXd& x, const MatrixXi& subsets, int& nkeep) {
	// initializations
	int n = x.size(), h = subsets.rows();
	VectorXi orderX = order(x), keep(nkeep);
	keep(0) = orderX(0);	// subset with lowest value of objective function
	// loop over subsets until nkeep unique subsets are found
	int j = 1, k = 1;
	while((k < nkeep) && (j < n)) {
		if(x(keep(k-1)) < x(orderX(j))) {
			// keep the current subset if its value for objective function is
			// larger than for the previous subset
			keep(k) = orderX(j);
			k++;
		} else {
			// equal value for objective function as for the previous subset
			// sort indices of previous and current subsets and check if they
			// are equal
			VectorXi previous = subsets.col(keep(k-1));
			sort(&previous.coeffRef(0), &previous.coeffRef(0)+h);
			VectorXi current = subsets.col(orderX(j));
			sort(&current.coeffRef(0), &current.coeffRef(0)+h);
			// loop over indices to see if they are all equal
			bool equal = true;
			int i = 0;
			while(equal && (i < h)) {
				equal = (previous(i) == current(i));
				i++;
			}
			// keep the current subset if it is not equal to the previous one
			if(!equal) {
				keep(k) = orderX(j);
				k++;
			}
		}
		j++;
	}
	// adjust nkeep if there are not enough unique subsets
	if(k < nkeep) {
		keep.conservativeResize(k);
		nkeep = k;
	}
	// return indices of best subsets
	return keep;
}

// compute the mean of a subset of the data
double subsetMean(const VectorXd& x, const VectorXi& subset) {
	const int h = subset.size();
	double mean = 0;
	for(int i = 0; i < h; i++) {
		mean += x(subset(i));
	}
	mean /= h;
	return mean;
}

// compute scale estimate on h smallest observations without correction factor
// (center estimate is passed as parameter)
double partialScale(const VectorXd& x, const double& center, const int& h) {
	// initialize STL vector for sorting
	const int n = x.size();
	vector<double> squares(n);
	for(int i = 0; i < n; i++) {
		squares[i] = pow(x(i)-center, 2);	// squared centered values
	}
	// call STL's nth_element()
	nth_element(squares.begin(), squares.begin()+h, squares.end());
	// compute scale estimate of h smallest observations
	double sumOfSquares = 0;
	for(int i = 0; i < h; i++) {
		sumOfSquares += squares[i];
	}
	return sqrt(sumOfSquares / double(h));
}


// sparse least trimmed squares
// Eigen library is used for linear algebra
// ATTENTION: intercept, coefficients, residuals, value of objective function,
//            residual center estimate and residual scale estimate are returned
//            through corresponding parameters
// x .............. predictor matrix
// y .............. response
// lambda ......... penalty parameter
// subsets ........ matrix of initial subsets
// useIntercept ... logical indicating whether intercept should be included
// ncstep ......... number of initial C-steps
// nkeep .......... number of subsets to keep after initial C-steps
// tol ............ numerical tolerance for convergence
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
// intercept ...... intercept is returned through this parameter
// coefficients ... coefficients are returned through this parameter
// residuals ...... residuals are returned through this parameter
// crit ........... value of objective function is returned through this parameter
// center ......... residual center estimate is returned through this parameter
// scale .......... residual scale estimate is returned through this parameter
VectorXi fastSparseLTS(const MatrixXd& x, const VectorXd& y,
		const double& lambda, Map<MatrixXi>& subsets, const bool& useIntercept,
		const int& ncstep, int& nkeep, const double& tol, const double& eps,
		const bool& useGram, double& intercept, VectorXd& coefficients,
		VectorXd& residuals, double& crit, double& center, double& scale) {
	// initializations
	const int n = x.rows(), p = x.cols();
	const int h = subsets.rows(), nsamp = subsets.cols();
	// compute lasso fits for initial subsets
	MatrixXd coefMat(p, nsamp);
	VectorXd intercepts;
	MatrixXd residMat(n, nsamp);
	VectorXd crits = VectorXd::Zero(nsamp);	// values of objective function
	if(useIntercept) {
		intercepts.resize(nsamp);
	}
	for(int k = 0; k < nsamp; k++) {
		// compute lasso fits
		coefMat.col(k) = fastLasso(x, y, lambda, subsets.col(k), useIntercept,
				eps, useGram, intercept, residuals, crit);
		residMat.col(k) = residuals;
		crits(k) = crit;
	}
	// perform initial c-steps on all subsets
	Array<bool, Dynamic, 1> continueCSteps(nsamp);
	for(int k = 0; k < nsamp; k++) {
		continueCSteps(k) = true;
	}
	for(int k = 0; k < nsamp; k++) {
		int i = 0;
		VectorXi subset = subsets.col(k);
		residuals = residMat.col(k);
		crit = crits(k);
		double previousCrit;
		while(continueCSteps(k) && (i < ncstep)) {
			previousCrit = crit;
			coefMat.col(k) = fastCStep(x, y, lambda, residuals, h, subset,
					useIntercept, eps, useGram, intercept, crit);
			continueCSteps(k) = ((previousCrit - crit) > tol);
			i++;
		}
		subsets.col(k) = subset;
		if(useIntercept) {
			intercepts(k) = intercept;
		}
		residMat.col(k) = residuals;
		crits(k) = crit;
	}

	// perform additional c-steps on best subsets until convergence and
	// keep track of optimal subset
	const VectorXi keep = bestSubsets(crits, subsets, nkeep);
	int which = keep(0);	// initial optimal subset
	for(int k = 0; k < nkeep; k++) {
		int keepK = keep(k);
		VectorXi subset = subsets.col(keepK);
		if(useIntercept) {
			intercept = intercepts(keepK);
		}
		residuals = residMat.col(keepK);
		crit = crits(keepK);
		double previousCrit;
		int i = ncstep;
		while(continueCSteps(keepK)) {
			previousCrit = crit;
			coefMat.col(keepK) = fastCStep(x, y, lambda, residuals, h, subset,
					useIntercept, eps, useGram, intercept, crit);
			continueCSteps(keepK) = ((previousCrit - crit) > tol);
			i++;
		}
		subsets.col(keepK) = subset;
		if(useIntercept) {
			intercepts(keepK) = intercept;
		}
		residMat.col(keepK) = residuals;
		crits(keepK) = crit;
		if(crit < crits(which)) {
			which = keepK;	// update optimal subset
		}
	}

	// return results
	VectorXi best = subsets.col(which);
	if(useIntercept) {
		intercept = intercepts(which);
	}
	coefficients = coefMat.col(which);
	residuals = residMat.col(which);
	crit = crits(which);
	center = subsetMean(residuals, best);  		// residual center
	scale = partialScale(residuals, center, h);	// uncorrected residual scale
	return best;
}

// R interface to fastSparseLTS()
// initial subsets are constructed in R and passed down to C++
SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subsets,
		SEXP R_intercept, SEXP R_ncstep, SEXP R_nkeep, SEXP R_tol, SEXP R_eps,
		SEXP R_useGram) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	double lambda = as<double>(R_lambda);
	IntegerMatrix Rcpp_subsets(R_subsets);	// matrix of initial subsets
	const int h = Rcpp_subsets.nrow(), nsamp = Rcpp_subsets.ncol();
	Map<MatrixXi> subsets(Rcpp_subsets.begin(), h, nsamp);	// reuse memory
	bool useIntercept = as<bool>(R_intercept);
	int ncstep = as<int>(R_ncstep);
	int nkeep = as<int>(R_nkeep);
	double tol = as<double>(R_tol);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function
	double intercept, crit, center, scale;
	VectorXd coefficients, residuals;
	VectorXi best = fastSparseLTS(x, y, lambda, subsets, useIntercept, ncstep,
			nkeep, tol, eps, useGram, intercept, coefficients, residuals,
			crit, center, scale);
	// currently only regression coefficients are returned
	NumericVector R_coefficients = wrap(coefficients);
	if(useIntercept) {
		R_coefficients.push_front(intercept);	// prepend intercept
	}
	return List::create(
			Named("best") = best,
			Named("coefficients") = R_coefficients,
			Named("residuals") = residuals,
			Named("crit") = crit,
			Named("center") = center,
			Named("scale") = scale
			);
}
