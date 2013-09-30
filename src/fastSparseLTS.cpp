/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#include "fastSparseLTS.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// *****************
// class definitions
// *****************

class Subset {
public:
	// information on the subset
	uvec indices;
	double intercept;
	vec coefficients;
	vec residuals;
	double crit;
	bool continueCSteps;

	// constructors
	Subset();
	Subset(const uword&, const uword&, const uword&);
	Subset(const uvec&);
	// compute lasso solution and residuals
	void lasso(const mat&, const vec&, const double&, const bool&,
			const bool&, const double&, const bool&);
	// perform C-Step
	void cStep(const mat&, const vec&, const double&, const bool&,
			const bool&, const double&, const double&, const bool&);
//	// overloaded < (is less) operator for sorting and ordering
//	bool operator< (const Subset&);
};

// constructors
inline Subset::Subset() {
	crit = R_PosInf;
	continueCSteps = true;
}
inline Subset::Subset(const uword&n, const uword&p, const uword& h) {
	indices = uvec(h);
	coefficients = vec(p);
	residuals = vec(n);
	crit = R_PosInf;
	continueCSteps = true;
}
inline Subset::Subset(const uvec& initial) {
	const uword h = initial.size();
	indices = uvec(h);
	for(uword i = 0; i < h; i++) {
		indices(i) = initial(i);	// data are copied
	}
	crit = R_PosInf;
	continueCSteps = true;
}

// compute lasso solution, residuals and value of objective function
void Subset::lasso(const mat& x, const vec& y, const double& lambda,
		const bool& normalize, const bool& useIntercept, const double& eps, 
    const bool& useGram) {
	// call standalone function
	fastLasso(x, y, lambda, true, indices, normalize, useIntercept, eps, 
      useGram, true, intercept, coefficients, residuals, crit);
}

// perform C-Step
void Subset::cStep(const mat& x, const vec& y, const double& lambda,
		const bool& normalize, const bool& useIntercept, const double& tol, 
    const double& eps, const bool& useGram) {
	// update subset
	const uword h = indices.size();
	indices = findSmallest(abs(residuals), h);
	// compute lasso solution for new subset
	// (also updated residuals and value of objective function)
	double previousCrit = crit;
	lasso(x, y, lambda, normalize, useIntercept, eps, useGram);
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

// compare two objects with == (is equal) operator
bool subsetIsEqual(const Subset& left, const Subset& right) {
	bool equal = (left.crit == right.crit);
	if(equal) {
		// values of the objective function are equal
		// check if indices are equal too
		uvec leftIndices = sort(left.indices);
		uvec rightIndices = sort(right.indices);
		// loop over sorted indices to see if they are all equal
		uword i = 0, h = leftIndices.size();
		while(equal && (i < h)) {
			equal = (leftIndices(i) == rightIndices(i));
			i++;
		}
	}
	return equal;
}


// R interface for object-based lasso (for testing purposes)
SEXP R_testLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
		SEXP R_normalize, SEXP R_intercept, SEXP R_eps, SEXP R_useGram) {
    // data initializations
	NumericMatrix Rcpp_x(R_x);						// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);		// reuse memory
	NumericVector Rcpp_y(R_y);			  // response
	vec y(Rcpp_y.begin(), n, false);	// reuse memory
	double lambda = as<double>(R_lambda);
	IntegerVector Rcpp_initial(R_initial);
	const int h = Rcpp_initial.size();
	uvec initial(h);
	for(int i = 0; i < h; i++) {
		initial(i) = Rcpp_initial[i] - 1;	// copy data
	}
  bool normalize = as<bool>(R_normalize);
	bool useIntercept = as<bool>(R_intercept);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// initialize object for subset and call native C++ function
	Subset subset(initial);
	subset.lasso(x, y, lambda, normalize, useIntercept, eps, useGram);
	// return results as list
	vec coefficients = subset.coefficients;
	if(useIntercept) {
		// prepend intercept
		coefficients.insert_rows(0, 1, false);
		coefficients(0) = subset.intercept;
	}
	return List::create(
			Named("indices") = subset.indices + 1,
			Named("coefficients") = coefficients,
			Named("residuals") = subset.residuals,
			Named("crit") = subset.crit,
			Named("continueCSteps") = subset.continueCSteps
			);
}

// R interface for object-based C-Step (for testing purposes)
SEXP R_testCStep(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subset,
		SEXP R_normalize, SEXP R_intercept, SEXP R_tol, SEXP R_eps, 
    SEXP R_useGram) {
    // data initializations
	NumericMatrix Rcpp_x(R_x);						// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);		// reuse memory
	NumericVector Rcpp_y(R_y);			  // response
	vec y(Rcpp_y.begin(), n, false);	// reuse memory
	double lambda = as<double>(R_lambda);
	List Rcpp_subset(R_subset);
  bool normalize = as<bool>(R_normalize);
	bool useIntercept = as<bool>(R_intercept);
	double tol = as<double>(R_tol);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// initialize object for subset
	IntegerVector Rcpp_indices = Rcpp_subset["indices"];
	const int h = Rcpp_indices.size();
	NumericVector Rcpp_coefficients = Rcpp_subset["coefficients"];
	NumericVector Rcpp_residuals = Rcpp_subset["residuals"];
	NumericVector Rcpp_crit = Rcpp_subset["crit"];
	Subset subset(n, p, h);
	for(int i = 0; i < h; i++) {
		subset.indices(i) = Rcpp_indices[i] - 1;
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
	subset.cStep(x, y, lambda, normalize, useIntercept, tol, eps, useGram);
	// return results as list
	vec coefficients = subset.coefficients;
	if(useIntercept) {
		// prepend intercept
		coefficients.insert_rows(0, 1, false);
		coefficients(0) = subset.intercept;
	}
	return List::create(
			Named("indices") = subset.indices + 1,
			Named("coefficients") = coefficients,
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
	// k counts the number of unique subsets found
	int k = 1, n = subsets.size();
	while((k < nkeep) && (k < n)) {
		if(subsetIsEqual(subsets[k-1], subsets[k])) {
			subsets.erase(subsets.begin()+k);
			n--;
		} else {
			k++;
		}
	}
	// adjust nkeep if there are not enough unique subsets
	if(k < nkeep) {
		nkeep = k;
	}
	// resize vector to keep the best subsets
	subsets.resize(nkeep);
}

// R interface for keeping best subsets (for testing purposes)
SEXP R_testKeepBest(SEXP R_subsetMat, SEXP R_crits, SEXP R_nkeep) {
	// data initializations
	IntegerMatrix Rcpp_subsetMat(R_subsetMat);		// subset matrix
	const int h = Rcpp_subsetMat.nrow(), nsamp = Rcpp_subsetMat.ncol();
	umat subsetMat(h, nsamp);
	for(int j = 0; j < nsamp; j++) {
		for(int i = 0; i < h; i++) {
			subsetMat(i,j) = Rcpp_subsetMat(i,j); // copy data
		}
	}
	NumericVector Rcpp_crits(R_crits);				    // values
	vec crits(Rcpp_crits.begin(), nsamp, false);	// reuse memory
	int nkeep = as<int>(R_nkeep);
	// call native C++ function
	vector<Subset> subsets(nsamp);
	for(int k = 0; k < nsamp; k++) {
		Subset subset(subsetMat.unsafe_col(k));
		subset.crit = crits(k);
		subsets[k] = subset;
	}
	keepBest(subsets, nkeep);
	// return results as list
	umat subsetMatOut(h, nkeep);
	vec critsOut(nkeep);
	for(int k = 0; k < nkeep; k++) {
		Subset subset = subsets[k];
		subsetMatOut.col(k) = subset.indices;
		critsOut(k) = subset.crit;
	}
	return List::create(
			Named("subsetMat") = subsetMatOut,
			Named("crits") = critsOut,
			Named("nkeep") = nkeep
			);
}

// compute the mean of a subset of the data
double subsetMean(const vec& x, const uvec& indices) {
//	return mean(x.elem(subset));
	const uword h = indices.size();
	double mean = 0;
	for(uword i = 0; i < h; i++) {
		mean += x(indices(i));
	}
	mean /= h;
	return mean;
}

// compute scale estimate on h smallest observations without correction factor
// (center estimate is passed as parameter)
double partialScale(const vec& x, const double& center, const int& h) {
	// initialize STL vector for sorting
	const int n = x.n_elem;
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
// Armadillo library is used for linear algebra
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
// ncores ......... number of processor cores for parallel computing
// center ......... residual center estimate is returned through this parameter
// scale .......... residual scale estimate is returned through this parameter
Subset fastSparseLTS(const mat& x, const vec& y, const double& lambda,
		const umat& initial, const bool& normalize, const bool& useIntercept, 
    const int& ncstep, int& nkeep, const double& tol, const double& eps, 
    const bool& useGram, int& ncores, double& center, double& scale) {
	// initializations
	const int h = initial.n_rows, nsamp = initial.n_cols;

	// block for parallel computing
	vector<Subset> subsets(nsamp);
	#pragma omp parallel num_threads(ncores)
	{
		// define STL vector of initial subsets, compute lasso fits and
		// perform initial c-steps
		#pragma omp for schedule(dynamic)
		for(int k = 0; k < nsamp; k++) {
			Subset subsetK(initial.unsafe_col(k));
			// compute lasso fit on initial subset
			subsetK.lasso(x, y, lambda, normalize, useIntercept, eps, useGram);
			// perform initial c-steps
			int i = 0;
			while(subsetK.continueCSteps && (i < ncstep)) {
				subsetK.cStep(x, y, lambda, normalize, useIntercept, tol, eps, useGram);
				i++;
			}
			subsets[k] = subsetK;
		}

		// keep best subsets (let only one process do the work)
		#pragma omp single
		if(nkeep < nsamp) {
			keepBest(subsets, nkeep);
		}

		// perform additional c-steps on best subsets until convergence
		#pragma omp for schedule(dynamic)
		for(int k = 0; k < nkeep; k++) {
			Subset subsetK = subsets[k];
			int i = 0;
			while(subsetK.continueCSteps) {
				subsetK.cStep(x, y, lambda, normalize, useIntercept, tol, eps, useGram);
				i++;
			}
			subsets[k] = subsetK;
		}
	}

	// find optimal subset
	int which = 0;					      // initial optimal subset
	double minCrit = R_PosInf;		// initial minimum of objective function
	for(int k = 0; k < nkeep; k++) {
		Subset subsetK = subsets[k];
		if(subsetK.crit < minCrit) {
			which = k;				      // update optimal subset
			minCrit = subsetK.crit;	// update minimum of objective function
		}
	}

	// return results
	Subset best = subsets[which];
	center = subsetMean(best.residuals, best.indices);  // residual center
	scale = partialScale(best.residuals, center, h);	  // uncorrected residual scale
	return best;
}

// R interface to fastSparseLTS()
// initial subsets are constructed in R and passed down to C++
SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
		SEXP R_normalize, SEXP R_intercept, SEXP R_ncstep, SEXP R_nkeep, 
    SEXP R_tol, SEXP R_eps, SEXP R_useGram, SEXP R_ncores) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);						// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);		// reuse memory
	NumericVector Rcpp_y(R_y);			  // response
	vec y(Rcpp_y.begin(), n, false);	// reuse memory
	double lambda = as<double>(R_lambda);
	IntegerMatrix Rcpp_initial(R_initial);	// matrix of initial subsets
	const int h = Rcpp_initial.nrow(), nsamp = Rcpp_initial.ncol();
	umat initial(h, nsamp);
	for(int j = 0; j < nsamp; j++) {
		for(int i = 0; i < h; i++) {
      // can't use the same memory-saving conversion for integer matrices
			initial(i,j) = Rcpp_initial(i,j) - 1;
		}
	}
  bool normalize = as<bool>(R_normalize);
  bool useIntercept = as<bool>(R_intercept);
	int ncstep = as<int>(R_ncstep);
	int nkeep = as<int>(R_nkeep);
	double tol = as<double>(R_tol);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	int ncores = as<int>(R_ncores);
	// call native C++ function
	double center, scale;
	Subset best = fastSparseLTS(x, y, lambda, initial, normalize, useIntercept,
			ncstep, nkeep, tol, eps, useGram, ncores, center, scale);
	// return results as list
	vec coefficients = best.coefficients;
	if(useIntercept) {
		// prepend intercept
		coefficients.insert_rows(0, 1, false);
		coefficients(0) = best.intercept;
	}
	return List::create(
			Named("best") = best.indices + 1,
			Named("coefficients") = coefficients,
			Named("residuals") = best.residuals,
			Named("objective") = best.crit,
			Named("center") = center,
			Named("scale") = scale
			);
}
