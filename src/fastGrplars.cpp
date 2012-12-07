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
// x ... predictor matrix
// y ... response variable
vec fastLm(const mat& x, const vec& y) {
	return solve(x, y);
}

// compute fitted values
// x ...... predictor matrix
// beta ... regression coefficients
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

//// find possible step sizes for groupwise LARS by solving quadratic equation
//vec computeStepSizes(const double& r, const double& a, const vec& corY,
//		const vec& corU, const vec& tau) {
//	// initializations
//	Environment robustHD("package:robustHD");
//	Function findStepSizes = robustHD["findStepSizes"];
//	NumericVector Rcpp_corY = wrap(corY.memptr(), corY.memptr() + corY.n_elem);
//	NumericVector Rcpp_corU = wrap(corU.memptr(), corU.memptr() + corU.n_elem);
//	NumericVector Rcpp_tau = wrap(tau.memptr(), tau.memptr() + tau.n_elem);
//	// call R function and convert result
//	NumericVector Rcpp_gammas = findStepSizes(r, a, corY, corU, tau);
//	vec gammas(Rcpp_gammas.begin(), Rcpp_gammas.size(), false);	// reuse memory
//	return gammas;
//}

// find smallest nonnegative real solution of a quadratic equation
double findSolution(const double& a, const double&b, const double& c) {
  // compute the discriminant and initialize the solution
  double discriminant = b*b - 4*a*c, solution = -b;
  if(discriminant > 0) {
    discriminant = sqrt(discriminant);
    vec solutions(2);
    solutions(0) = solution + discriminant;
    solutions(1) = solution - discriminant;
    solutions /= 2*a;  // divide by denominator
    // keep only nonnegative solutions
    const uvec keep = find(solutions >= 0);
    solutions = solutions.elem(keep);
    // check whether there is a nonnegative solution
    const uword n = solutions.n_elem;
    if(n == 0) {
      solution = R_PosInf;
    } else if(n == 1) {
      solution = solutions(0);
    } else {
      solution = solutions.min(); // take the smallest one
    }
  } else {
    solution /= 2*a;  // divide by denominator
    // check whether the solution is nonnegative
    if(solution < 0) solution = R_PosInf;
  }
  // return the smallest nonnegative solution
  return solution;
}

// find possible step sizes for groupwise LARS by solving quadratic equation
vec computeStepSizes(const double& r, const double& a, const vec& corY,
  	const vec& corU, const vec& tau) {
	// initializations
  const uword n = corY.n_elem;
  vec gammas(n);
  // compute step size for each predictor group
  for(uword j = 0; j < n; j++) {
    gammas(j) = findSolution(a*a - tau(j)*tau(j), 2 * (corY(j)*corU(j) - r*a), 
        r*r - corY(j)*corY(j));
  }
  // return step sizes
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
// ncores ... number of processor cores for parallel computing
// parallel computing is only used for expensive computations with all
// inactive predictors, otherwise there is no speedup due to overhead
uvec fastGrplars(const mat& x, const vec& y, const uword& sMax,
		const vector<uvec>& assign, int& ncores) {
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
		if(p(j) != p(0)) {
			adjust = true;
			break;
		}
	}

	// STEP 1: find first ranked predictor group
	// compute the correlation of the fitted values from each predictor group
	// with the response
	mat yHat(n, m);
	vec corY(m);
	#pragma omp parallel for num_threads(ncores) schedule(dynamic)
	for(uword j = 0; j < m; j++) {
		mat xj = x.cols(assign[j]);
		vec beta = fastLm(xj, y);
		yHat.col(j) = fitted(xj, beta);
		corY(j) = stddev(yHat.unsafe_col(j));
	}
	// if necessary, compute the denominators for adjustment w.r.t. the number
	// of variables in the group, and adjust the correlations of the fitted
	// values with the response
	vec adjustment;
	if(adjust) {
		adjustment.set_size(m);
		for(uword j = 0; j < m; j++) {
			adjustment(j) = sqrt((double)(p(j)));
			corY(j) /= adjustment(j);
		}
	}
	// find predictor group with maximum correlation
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
	// initialize correlation matrix of the fitted values from the active
    // predictor groups
	mat R = eye<mat>(1, 1);
	// initialize correlation of the fitted values from the active predictors
	// with the equiangular vector
	double a = 1.0;
	if(adjust) {
		a /= adjustment(whichMax);	// adjust for unequal group sizes
	}
	// other quantities related to the equiangular vector
	vec w = ones<vec>(1), u(n), corU(m-1), tau(m-1);
	// keep track of the scale of the current response for update formulas
	double sigma = 1.0;
	// start iterative computations
	for(uword k = 1; k < sMax; k++) {
		// standardize the fitted values for the new active predictor group
		standardize(yHat, active(k-1));
		// compute the equiangular vector
		if(k == 1) {
			u = yHat.col(whichMax);	// fitted values from the first active group
		} else {
			// expand correlation matrix of the fitted values from the active
			// predictor groups
			R.resize(k, k);
			R(k-1, k-1) = 1;
			vec yk = yHat.unsafe_col(active(k-1));
			for(uword j = 0; j < k-1; j++) {
				R(j, k-1) = corPearson(yk, yHat.unsafe_col(active(j)));
				R(k-1, j) = R(j, k-1);
			}
			// other computations according to algorithm
			mat invR = solve(R, eye<mat>(k, k));
			vec q(k);
			if(adjust) {
				q = adjustment.elem(active);
			} else {
				q.ones();
			}
			// compute the correlation of the fitted values from the active
			// predictor groups with the equiangular vector (adjustment for
			// unequal group size is considered via vector 'q')
			a = 1 / sqrt(as_scalar(conv_to<rowvec>::from(q) * invR * q));
			// compute the equiangular vector
			w = a * (invR * q);
			u = yHat.cols(active) * w;
		}
		// compute the fitted values of the equiangular vector for each
		// inactive predictor group, as well as the correlations involving the
		// inactive predictor groups and the equiangular vector
		mat uHat(n, inactive.n_elem);
		#pragma omp parallel for num_threads(ncores) schedule(dynamic)
		for(uword j = 0; j < inactive.n_elem; j++) {
			// compute the fitted values
			mat xj = x.cols(assign[inactive(j)]);
			vec beta = fastLm(xj, u);
			uHat.col(j) = fitted(xj, beta);
			// compute the correlations
			corU(j) = corPearson(yHat.unsafe_col(inactive(j)), u);
			tau(j) = stddev(uHat.unsafe_col(j));
		}
		if(adjust) {
			// adjustment for unequal group size
			for(uword j = 0; j < inactive.n_elem; j++) {
				double tmp = adjustment(inactive(j));
				corU(j) /= tmp;
				tau(j) /= tmp;
			}
		}
		// compute the step size by solving the quadratic equation
		vec gammas = computeStepSizes(r(k-1), a, corY, corU, tau);
		uword whichMin;
		double gamma = gammas.min(whichMin);
        // the following computations are not necessary in the last iteration
		if(k < (sMax-1)) {
    		// update the scale of the current response
            sigma = sqrt(1 - 2 * gamma * r(k-1)/a + pow(gamma, 2));
           	// update the fitted values from the new or not yet sequenced
            // predictor groups
			for(uword j = 0; j < inactive.n_elem; j++) {
        		vec yj = yHat.unsafe_col(inactive(j));
        		yj -= gamma * uHat.unsafe_col(j);
        		yj /= sigma;
        	}
			// update correlations (adjustment for unequal group size is taken
        	// care of by update formula)
			r.insert_rows(k, 1, false);	// do not initialize new memory
			r(k) = (r(k-1) - gamma * a) / sigma;
 			corY.shed_row(whichMin);
			corU.shed_row(whichMin);
			tau.shed_row(whichMin);
			// this is not alias safe:
//			corY = sqrt(pow(corY, 2) - 2 * gamma * corU * corY +
//					pow(gamma, 2) * pow(tau, 2)) / sigma;
			// this works:
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
SEXP R_fastGrplars(SEXP R_x, SEXP R_y, SEXP R_sMax,
		SEXP R_assign, SEXP R_ncores) {
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
	int ncores = as<int>(R_ncores);
	// call native C++ function
	uvec active = fastGrplars(x, y, sMax, assign, ncores);
	return wrap(active.memptr(), active.memptr() + active.n_elem);
}
