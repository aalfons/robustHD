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

// standardize data in place using mean and standard deviation
// x ........ data matrix
// select ... indices of columns to be standardized
void standardize(const mat& x, const uword& select) {
	const uword n = x.n_rows;
	// with unsafe_col(), standardization is done in place
	vec xj = x.unsafe_col(select);
	double center = mean(xj);							// compute mean
	xj -= center;										// sweep out mean
	double scale = norm(xj, 2) / sqrt((double)(n-1));	// compute SD
	xj /= scale;										// sweep out SD
}

// find possible step sizes for groupwise LARS by solving quadratic equation
vec computeStepSizes(const double& r, const double& a, const vec& corY,
		const vec& corU, const vec& tau) {
	// initializations
	Environment robustHD("package:robustHD");
	Function findStepSizes = robustHD["findStepSizes"];
	NumericVector Rcpp_corY = wrap(corY.memptr(), corY.memptr() + corY.n_elem);
	NumericVector Rcpp_corU = wrap(corU.memptr(), corU.memptr() + corU.n_elem);
	NumericVector Rcpp_tau = wrap(tau.memptr(), tau.memptr() + tau.n_elem);
	// call R function and convert result
	NumericVector Rcpp_gammas = findStepSizes(r, a, corY, corU, tau);
	vec gammas(Rcpp_gammas.begin(), Rcpp_gammas.size(), false);	// reuse memory
	return gammas;
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
	// compute correlation with response for each predictor group
	mat yHat(n, m);
	vec corY(m);
	for(uword j = 0; j < m; j++) {
		mat xj = x.cols(assign[j]);
		vec beta = fastLm(xj, y);
		yHat.col(j) = fitted(xj, beta);
		corY(j) = stddev(yHat.unsafe_col(j));
	}
	vec sqrtP;
	if(adjust) {
		sqrtP.set_size(m);
		// compute denominators for adjustment w.r.t. the number of variables
		// in the group and adjust correlation with response
		for(uword j = 0; j < m; j++) {
			sqrtP(j) = sqrt((double)(p(j)));
			corY(j) /= sqrtP(j);
		}
	}
	// find predictor with maximum correlation
	vec r(1);
	uword whichMax;
	r(0) = corY.max(whichMax);	// find index of maximum correlation
	// initialize active set
	uvec active(1);
	active(0) = whichMax;
	// initialize inactive set
	uvec inactive = seqLen(m);
	inactive.shed_row(whichMax);
    corY.shed_row(whichMax);

	// STEP 2: update active set
	// correlation matrix of active predictor groups
	mat R = eye<mat>(1, 1);
	// quantities related to the equiangular vector
	double a = 1.0;
	if(adjust) {
		// adjust correlation of active variables equiangular vector
		a /= sqrtP(whichMax);
	}
	vec w = ones<vec>(1), u(n), corU(m-1), tau(m-1);
	// keep track of scale for update formulas
	double sigma = 1.0;
	// start iterative computations
	for(uword k = 1; k < sMax; k++) {
		// standardize fitted values for new active predictor group
		standardize(yHat, active(k-1));
		// compute the equiangular vector
		if(k == 1) {
			u = yHat.col(whichMax);
		} else {
			// compute correlations between fitted values for active groups
			R.resize(k, k);
			R(k-1, k-1) = 1;
			vec yk = yHat.unsafe_col(active(k-1));
			for(uword j = 0; j < k-1; j++) {
				R(j, k-1) = corPearson(yk, yHat.unsafe_col(active(j)));
				R(k-1, j) = R(j, k-1);
			}
			// other computations according to algorithm
			mat invR = solve(R, eye<mat>(k, k));
			vec q;
			if(adjust) {
				q = sqrtP.elem(active);
			} else {
				q.ones(k);
			}
			// compute correlation of active predictor groups with equiangular
			// vector (adjustment for unequal group size is considered via
			// vector 'q')
			a = 1 / sqrt(as_scalar(trans(q) * invR * q));
			// compute equiangular vector
			w = a * (invR * q);
			u = yHat.cols(active) * w;
		}
		// compute the fitted values of the equiangular vector for each
		// inactive predictor group
		mat uHat(n, inactive.n_elem);
		for(uword j = 0; j < inactive.n_elem; j++) {
			mat xj = x.cols(assign[inactive(j)]);
			vec beta = fastLm(xj, u);
			uHat.col(j) = fitted(xj, beta);
		}
		// compute correlations involving the inactive predictor groups and
		// the equiangular vector
		for(uword j = 0; j < inactive.n_elem; j++) {
			corU(j) = corPearson(yHat.unsafe_col(inactive(j)), u);
			tau(j) = stddev(uHat.unsafe_col(j));
		}
		if(adjust) {
			// adjustment for unequal group size
			for(uword j = 0; j < inactive.n_elem; j++) {
				double tmp = sqrtP(inactive(j));
				corU(j) /= tmp;
				tau(j) /= tmp;
			}
		}
		// compute step size by solving quadratic equation
		vec gammas = computeStepSizes(r(k-1), a, corY, corU, tau);
		uword whichMin;
		double gamma = gammas.min(whichMin);
        // the following computations are not necessary in the last iteration
		if(k < (sMax-1)) {
    		// update scale of response
            sigma = sqrt(1 - 2 * gamma * r(k-1)/a + pow(gamma, 2));
           	// update fitted values for new or not yet sequenced groups
        	for(uword j = 0; j < inactive.size(); j++) {
        		vec tmp = yHat.unsafe_col(inactive(j));
        		tmp -= gamma * uHat.unsafe_col(j);
        		tmp /= sigma;
        	}
			// update correlations (adjustment for unequal group size is taken
        	// care of by update formula)
			r.insert_rows(k, 1, false);	// do not initialize new memory
			r(k) = (r(k-1) - gamma * a) / sigma;
 			corY.shed_row(whichMin);
			corU.shed_row(whichMin);
			tau.shed_row(whichMin);
			// this is not alias safe:
//
//			corY = sqrt(pow(corY, 2) - 2 * gamma * corU * corY +
//					pow(gamma, 2) * pow(tau, 2)) / sigma;
			// this works fine:
			for(uword j = 0; j < corY.n_elem; j++) {
				corY(j) = sqrt(pow(corY(j), 2) - 2*gamma*corU(j)*corY(j) +
						pow(gamma, 2) * pow(tau(j), 2)) / sigma;
			}
        }
        // update active set
		active.insert_rows(k, 1, false);	// do not initialize new memory
		active(k) = inactive(whichMin);
		// update inactive set
		inactive.shed_row(whichMin);
	}

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
