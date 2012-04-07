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
	VectorXi indices;
	double intercept;
	VectorXd coefficients;
	VectorXd residuals;
	double crit;
	bool continueCSteps;

	// constructors
	Subset();
	Subset(const int&, const int&, const int&);
	Subset(const VectorXi&);
	// compute lasso solution and residuals
	void lasso(const MatrixXd&, const VectorXd&, const double&, const bool&,
			const double&, const bool&);
	// compute objective function
	void objective(const double&);
	// perform C-Step
	void cStep(const MatrixXd&, const VectorXd&, const double&, const bool&,
			const double&, const double&, const bool&);
//	// overloaded < (is less) operator for sorting and ordering
//	bool operator< (const Subset&);
};

// constructors
inline Subset::Subset() {
	crit = R_PosInf;
	continueCSteps = true;
}
inline Subset::Subset(const int&n, const int&p, const int& h) {
	indices = VectorXi(h);
	coefficients = VectorXd(p);
	residuals = VectorXd(n);
	crit = R_PosInf;
	continueCSteps = true;
}
inline Subset::Subset(const VectorXi& initial) {
	const int h = initial.size();
	indices = VectorXi(h);
	for(int i = 0; i < h; i++) {
		indices(i) = initial(i);	// data are copied
	}
	crit = R_PosInf;
	continueCSteps = true;
}

// compute sparse LTS objective function
// (L1 penalized trimmed sum of squared residuals)
void Subset::objective(const double& lambda) {
	// compute sum of squared residuals for subset
	const int h = indices.size();
	crit = 0;
	for(int i = 0; i < h; i++) {
		crit += pow(residuals(indices(i)), 2);
	}
	// add L1 penalty on coefficients
	crit += h * lambda * coefficients.lpNorm<1>();
}

// compute lasso solution, residuals and value of objective function
void Subset::lasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useIntercept, const double& eps, const bool& useGram) {
	// compute coefficients
	coefficients = fastLasso(x, y, lambda, true, indices,
			useIntercept, eps, useGram, intercept);
	// compute residuals
	residuals.noalias() = y - x * coefficients;
	if(useIntercept) {
		const int n = residuals.size();
		for(int i = 0; i < n; i++) {
			residuals(i) -= intercept;
		}
	}
	// compute value of objective function
	objective(lambda);
}

// perform C-Step
void Subset::cStep(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useIntercept, const double& tol, const double& eps,
		const bool& useGram) {
	// update subset
	const int h = indices.size();
	indices = findSmallest(residuals.cwiseAbs(), h);
	// compute lasso solution for new subset
	// (also updated residuals and value of objective function)
	double previousCrit = crit;
	lasso(x, y, lambda, useIntercept, eps, useGram);
	// check for convergence
	continueCSteps = ((previousCrit - crit) > tol);
}

//// overloaded < (is less) operator for sorting and ordering
//bool Subset::operator< (const Subset& other) {
//      return (this->value < other.value);
//}

// compare two objects with < (is less) operator for sorting and ordering
bool subsetIsLess(const Subset& left, const Subset& right) {
	return left.crit < right.crit;
}


// R interface for object-based lasso (for testing purposes)
SEXP R_testLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
		SEXP R_intercept, SEXP R_eps, SEXP R_useGram) {
    // data initializations
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	double lambda = as<double>(R_lambda);
	IntegerVector Rcpp_initial(R_initial);
	const int h = Rcpp_initial.size();
	VectorXi initial(h);
	for(int i = 0; i < h; i++) {
		initial(i) = Rcpp_initial[i];	// do not reuse memory
	}
	bool useIntercept = as<bool>(R_intercept);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// initialize object for subset and call native C++ function
	Subset subset(initial);
	subset.lasso(x, y, lambda, useIntercept, eps, useGram);
	// return results as list
	NumericVector Rcpp_coefficients = wrap(subset.coefficients);
	if(useIntercept) {
		Rcpp_coefficients.push_front(subset.intercept);	// prepend intercept
	}
	return List::create(
			Named("indices") = subset.indices,
			Named("coefficients") = Rcpp_coefficients,
			Named("residuals") = subset.residuals,
			Named("crit") = subset.crit,
			Named("continueCSteps") = subset.continueCSteps
			);
}

// R interface for object-based C-Step (for testing purposes)
SEXP R_testCStep(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subset,
		SEXP R_intercept, SEXP R_tol, SEXP R_eps, SEXP R_useGram) {
    // data initializations
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	double lambda = as<double>(R_lambda);
	List Rcpp_subset(R_subset);
	bool useIntercept = as<bool>(R_intercept);
	double tol = as<double>(R_tol);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// initialize object for subset
	NumericVector Rcpp_indices = Rcpp_subset["indices"];
	const int h = Rcpp_indices.size();
	NumericVector Rcpp_coefficients = Rcpp_subset["coefficients"];
	NumericVector Rcpp_residuals = Rcpp_subset["residuals"];
	NumericVector Rcpp_crit = Rcpp_subset["crit"];
	Subset subset(n, p, h);
	for(int i = 0; i < h; i++) {
		subset.indices(i) = Rcpp_indices[i];
	}
	if(useIntercept) {
		subset.intercept = Rcpp_coefficients[0];
		Rcpp_coefficients.erase(0);
	}
	for(int j = 0; j < p; j++) {
		subset.coefficients(j) = Rcpp_coefficients[j];
	}
	for(int i = 0; i < n; i++) {
		subset.residuals(i) = Rcpp_residuals[i];
	}
	subset.crit = Rcpp_crit[0];
	// call native C++ function
	subset.cStep(x, y, lambda, useIntercept, tol, eps, useGram);
	// return results as list
	Rcpp_coefficients = wrap(subset.coefficients);
	if(useIntercept) {
		Rcpp_coefficients.push_front(subset.intercept);	// prepend intercept
	}
	return List::create(
			Named("indices") = subset.indices,
			Named("coefficients") = Rcpp_coefficients,
			Named("residuals") = subset.residuals,
			Named("crit") = subset.crit,
			Named("continueCSteps") = subset.continueCSteps
			);
}


// ****************************
// sparse least trimmed squares
// ****************************

// keep best subsets
void keepBest(vector<Subset>& subsets, int& nkeep) {
	// sort subsets
	sort(subsets.begin(), subsets.end(), subsetIsLess);
	// loop over subsets until nkeep unique subsets are found
	int k = 1, n = subsets.size();
	while((k < nkeep) && (k < n)) {
		if(subsetIsLess(subsets[k-1], subsets[k])) {
			// keep the current subset if its value for objective function is
			// larger than for the previous subset
			k++;
		} else {
			// equal value for objective function as for the previous subset
			// sort indices of previous and current subsets and check if they
			// are equal
			VectorXi previous = (subsets[k-1]).indices;
			VectorXi current = (subsets[k]).indices;
			int h = previous.size();
			sort(&previous.coeffRef(0), &previous.coeffRef(0)+h);
			sort(&current.coeffRef(0), &current.coeffRef(0)+h);
			// loop over indices to see if they are all equal
			bool equal = true;
			int i = 0;
			while(equal && (i < h)) {
				equal = (previous(i) == current(i));
				i++;
			}
			// remove the current subset if it is equal to the previous one,
			// otherwise keep it
			if(equal) {
				subsets.erase(subsets.begin()+k);
				n--;
			} else {
				k++;
			}
		}
	}
	// adjust nkeep if there are not enough unique subsets
	if(k < nkeep) {
		nkeep = k;
	}
	// resize vector to keep the best subsets
	subsets.resize(nkeep);
}

// compute the mean of a subset of the data
double subsetMean(const VectorXd& x, const VectorXi& indices) {
	const int h = indices.size();
	double mean = 0;
	for(int i = 0; i < h; i++) {
		mean += x(indices(i));
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
// initial ........ matrix of initial subsets
// useIntercept ... logical indicating whether intercept should be included
// ncstep ......... number of initial C-steps
// nkeep .......... number of subsets to keep after initial C-steps
// tol ............ numerical tolerance for convergence
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
// center ......... residual center estimate is returned through this parameter
// scale .......... residual scale estimate is returned through this parameter
Subset fastSparseLTS(const MatrixXd& x, const VectorXd& y,
		const double& lambda, const MatrixXi& initial, const bool& useIntercept,
		const int& ncstep, int& nkeep, const double& tol, const double& eps,
		const bool& useGram, double& center, double& scale) {
	// initializations
	const int n = x.rows(), p = x.cols();
	const int h = initial.rows(), nsamp = initial.cols();
	// define STL vector of initial subsets and compute lasso fits
	vector<Subset> subsets(nsamp);
	for(int k = 0; k < nsamp; k++) {
		Subset subsetK(initial.col(k));
		subsetK.lasso(x, y, lambda, useIntercept, eps, useGram);
		subsets[k] = subsetK;
	}
	// perform initial c-steps on all subsets
	for(int k = 0; k < nsamp; k++) {
		Subset subsetK = subsets[k];
		int i = 0;
		while(subsetK.continueCSteps && (i < ncstep)) {
			subsetK.cStep(x, y, lambda, useIntercept, tol, eps, useGram);
			i++;
		}
		subsets[k] = subsetK;
	}

	// perform additional c-steps on best subsets until convergence and
	// keep track of optimal subset
	if(nkeep < nsamp) {
		keepBest(subsets, nkeep);	// keep best subsets
	}
	int which = 0;					// initial optimal subset
	double minCrit = R_PosInf;		// initial minimum of objective function
	for(int k = 0; k < nkeep; k++) {
		Subset subsetK = subsets[k];
		int i = 0;
		while(subsetK.continueCSteps) {
			subsetK.cStep(x, y, lambda, useIntercept, tol, eps, useGram);
			i++;
		}
		if(subsetK.crit < minCrit) {
			which = k;	// update optimal subset
		}
		subsets[k] = subsetK;
	}

	// return results
	Subset best = subsets[which];
	center = subsetMean(best.residuals, best.indices);  // residual center
	scale = partialScale(best.residuals, center, h);	// uncorrected residual scale
	return best;
}

// R interface to fastSparseLTS()
// initial subsets are constructed in R and passed down to C++
SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
		SEXP R_intercept, SEXP R_ncstep, SEXP R_nkeep, SEXP R_tol, SEXP R_eps,
		SEXP R_useGram) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	double lambda = as<double>(R_lambda);
	IntegerMatrix Rcpp_initial(R_initial);	// matrix of initial subsets
	const int h = Rcpp_initial.nrow(), nsamp = Rcpp_initial.ncol();
	Map<MatrixXi> initial(Rcpp_initial.begin(), h, nsamp);	// reuse memory
	bool useIntercept = as<bool>(R_intercept);
	int ncstep = as<int>(R_ncstep);
	int nkeep = as<int>(R_nkeep);
	double tol = as<double>(R_tol);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function
	double center, scale;
	Subset best = fastSparseLTS(x, y, lambda, initial, useIntercept, ncstep,
			nkeep, tol, eps, useGram, center, scale);
	// return results as list
	NumericVector Rcpp_coefficients = wrap(best.coefficients);
	if(useIntercept) {
		Rcpp_coefficients.push_front(best.intercept);	// prepend intercept
	}
	return List::create(
			Named("best") = best.indices,
			Named("coefficients") = Rcpp_coefficients,
			Named("residuals") = best.residuals,
			Named("crit") = best.crit,
			Named("center") = center,
			Named("scale") = scale
			);
}
