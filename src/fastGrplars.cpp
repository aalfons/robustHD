/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include <R.h>
#include "fastGrplars.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// barebones version of least squares regression
vec fastLm(const mat& x, const vec& y) {
	return solve(x, y);
}

// compute fitted values
vec fitted(const mat& x, const vec& beta) {
	return x * beta;
}


// variable sequencing via robust least angle regression
// Armadillo library is used for linear algebra
// ATTENTION: the data are assumed to be standardized (this is done in R)
// x ........ predictor matrix
// y ........ response
// sMax ..... number of predictors to be sequenced
// assign ... each element is an integer vector giving the variables that
//            belong to the corresponding group
uvec fastGrplars(const mat& x, const vec& y, const uword& sMax,
		const vector<uvec>& assign) {
	// initializations
	const uword n = x.n_rows, m = assign.size();
	// determine number of variables in each predictor group
	uvec p(m);
	for(uword j = 0; j < m; j++) {
		p(j) = (assign[j]).n_elem;
	}
	// determine whether adjustment for different group sizes is necessary
	bool adjust = false;
	for(uword j = 1; j < m; j++) {
		if(p(j) != p(1)) {
			adjust = true;
			break;
		}
	}

	// STEP 1: find first ranked predictor group
	// compute R-squared for each predictor group
	mat yHat(n, m);
	vec Rsq(m);
	for(uword j = 0; j < m; j++) {
		mat xj = x.cols(assign[j]);
		vec beta = fastLm(xj, y);
		yHat.col(j) = fitted(xj, beta);
		Rsq(j) = var(yHat.unsafe_col(j));
	}
	// find predictor with maximum R-squared
	uword whichMax;
	if(adjust) {
		// adjusted R-squared w.r.t. the number of variables in the group
		vec adjustedRsq(m);
		for(uword j = 0; j < m; j++) {
			adjustedRsq(j) = Rsq(j) / (double)(p(j));
		}
		// find index of maximum adjusted R-squared
		adjustedRsq.max(whichMax);
	} else {
		Rsq.max(whichMax);	// find index of maximum R-squared
	}
	// initialize active set
	uvec active(1);
	active(0) = whichMax;
	// initialize inactive set
	uvec inactive = seqLen(m);
	inactive.shed_row(whichMax);

	// STEP 2: update active set

	// return active set
	return active;
}

// R interface to fastRlars()
SEXP R_fastGrplars(SEXP R_x, SEXP R_y, SEXP R_sMax, SEXP R_assign) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);						// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);				// reuse memory
	NumericVector Rcpp_y(R_y);			// response
	vec y(Rcpp_y.begin(), n, false);	// reuse memory
	uword sMax = as<uword>(R_sMax);
	// convert list giving group assignment of predictors to C++ data structure
	List Rcpp_assign(R_assign);
	const int m = Rcpp_assign.size();
	vector<uvec> assign(m);
	for(int j = 0; j < m; j++) {
		SEXP R_group = Rcpp_assign[j];
		IntegerVector Rcpp_group(R_group);
		const int pj = Rcpp_group.size();
		uvec group(pj);
		for(int k = 0; k < pj; k++) {
			group(k) = Rcpp_group[k] - 1;
		}
		assign[j] = group;	// indices of variables in j-th predictor group
	}
	// call native C++ function
	uvec active = fastGrplars(x, y, sMax, assign);
	return wrap(active.memptr(), active.memptr() + active.n_elem);
}
