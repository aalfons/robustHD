/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include <R.h>
#include "fastRlars.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// apply scale function to columns of a matrix
vec applyScaleFun(const mat& x, SEXP scaleFun) {
	// initializations
	Environment base("package:base");
	Function apply = base["apply"];
	NumericMatrix Rcpp_x = wrap(x);	// does this reuse memory?
	// call R function and convert result
	NumericVector Rcpp_scale = apply(Rcpp_x, 2, scaleFun);
	vec scale(Rcpp_scale.begin(), Rcpp_scale.size(), false);	// reuse memory
	return scale;
}


// variable sequencing via robust least angle regression
// Armadillo library is used for linear algebra
// ATTENTION: the data are assumed to be standardized (this is done in R)
// x .......... predictor matrix
// y .......... response
// sMax ....... number of predictors to be sequenced
// c .......... tuning constant for initial adjusted univariate winsorization
// prob ....... tuning parameter for bivariate winsorization
// tol ........ numeric tolerance for detecting singularity
// scaleFun ... R function to compute scale estimates
// ncores ..... number of processor cores for parallel computing
// parallel computing is only used for expensive computations with all
// inactive predictors, otherwise there is no speedup due to overhead
uvec fastRlars(const mat& x, const vec& y, const uword& sMax, const double& c,
		const double& prob, const double& tol, SEXP scaleFun, int& ncores) {
	// initializations
	const uword p = x.n_cols;

	// STEP 1: find first ranked predictor
	// compute correlations with response
	vec corY(p);
	#pragma omp parallel for num_threads(ncores) schedule(dynamic)
	for(uword j = 0; j < p; j++) {
		corY(j) = corHuberBi(x.unsafe_col(j), y, c, prob, tol);
	}
	vec absCorY = abs(corY);
	// find predictor with maximum absolute correlation
	vec r(1);
	ivec signs(1);
	uword whichMax;
	r(0) = absCorY.max(whichMax);		// absolute correlation of active predictor
	signs(0) = sign(corY(whichMax));	// sign of correlation of active predictor
	// initialize active set
	uvec active(1);
	active(0) = whichMax;
	// initialize inactive set
	uvec inactive = seqLen(p);
	inactive.shed_row(whichMax);
    corY.shed_row(whichMax);
    uword m = inactive.n_elem;

	// STEP 2: update active set
	// further initializations
	mat R = ones<mat>(p, sMax);	// correlation matrix of predictors with active set
	double a = 1.0;				// correlations of active variables with equiangular vector
    vec w = ones<vec>(1);		// coefficients of active variables for equiangular vector
	// start iterative computations
	for(uword k = 1; k < sMax; k++) {
		// compute correlations of inactive predictors with new active predictor
		vec xk = x.unsafe_col(active(k-1));
		#pragma omp parallel for num_threads(ncores) schedule(dynamic)
		for(uword j = 0; j < m; j++) {
			vec xj = x.unsafe_col(inactive(j));
			R(inactive(j), k-1) = corHuberBi(xj, xk, c, prob, tol);
		}
		for(uword j = 1; j < k; j++) {
			R(active(j-1), k-1) = R(active(k-1), j-1);
		}
		// compute correlations of active variables with equiangular vector
        if(k > 1) {
        	// correlation matrix taking signs into account
        	mat G(k, k);
        	for(uword j = 0; j < k; j++) {
        		for(uword i = 0; i < k; i++) {
        			G(i, j) = signs(i) * signs(j) * R(active(i), j);
        		}
        	}
        	// check if correlation matrix is positive definite
        	// compute eigenvalues and eigenvectors
        	vec eigVal;
        	mat eigVec;
        	eig_sym(eigVal, eigVec, G);
        	if(eigVal(0) < 0) {  // first eigenvalue is the smallest
        		// correction of correlation matrix for positive definiteness
        		vec lambda = square(applyScaleFun(x.cols(active) * eigVec, scaleFun));
        		G = eigVec * diagmat(lambda) * trans(eigVec);
        	}
            // compute quantities necessary for computing the step size
        	mat invG = solve(G, eye<mat>(k, k));
            a = 1 / sqrt(as_scalar(ones<rowvec>(k) * invG * ones<vec>(k)));
            w = a * (invG * ones<vec>(k));
        }
        // compute correlations of inactive predictors with equiangular vector
        vec corU(m);
        #pragma omp parallel for num_threads(ncores) schedule(dynamic)
        for(uword j = 0; j < m; j++) {
        	corU(j) = sum(signs % R(inactive(j), span(0, k-1)) % w);
        }
        // compute step size in equiangular direction
        vec gammaMinus = (r(k-1) - corY) / (a - corU);
        vec gammaPlus = (r(k-1) + corY) / (a + corU);
        for(uword j = 0; j < m; j++) {
        	if(gammaMinus(j) <= 0) gammaMinus(j) = R_PosInf;
        	if(gammaPlus(j) <= 0) gammaPlus(j) = R_PosInf;
        }
        uword whichMinus, whichPlus, whichMin;
        double minGammaMinus, minGammaPlus, gamma;
        minGammaMinus = gammaMinus.min(whichMinus);
        minGammaPlus = gammaPlus.min(whichPlus);
		signs.insert_rows(k, 1, false);	// do not initialize new memory
        if(minGammaMinus < minGammaPlus) {
        	whichMin = whichMinus;
        	gamma = minGammaMinus;
        	signs(k) = 1;
        } else {
        	whichMin = whichPlus;
        	gamma = minGammaPlus;
        	signs(k) = -1;
        }
        // update correlations
		r.insert_rows(k, 1, false);			// do not initialize new memory
        r(k) = r(k-1) - gamma * a;
		corY.shed_row(whichMin);
		corU.shed_row(whichMin);
		corY = corY - gamma * corU;
		// update active set
		active.insert_rows(k, 1, false);	// do not initialize new memory
		active(k) = inactive(whichMin);
		// update inactive set
		inactive.shed_row(whichMin);
		m--;	// decrement number of inactive variables
	}

	// return active set
	return active;
}

// R interface to fastRlars()
SEXP R_fastRlars(SEXP R_x, SEXP R_y, SEXP R_sMax, SEXP R_c, SEXP R_prob,
		SEXP R_tol, SEXP scaleFun, SEXP R_ncores) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);						// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);				// reuse memory
	NumericVector Rcpp_y(R_y);			// response
	vec y(Rcpp_y.begin(), n, false);	// reuse memory
	uword sMax = as<uword>(R_sMax);
	double c = as<double>(R_c);
	double prob = as<double>(R_prob);
	double tol = as<double>(R_tol);
	int ncores = as<int>(R_ncores);
	// call native C++ function
	uvec active = fastRlars(x, y, sMax, c, prob, tol, scaleFun, ncores);
	return wrap(active.memptr(), active.memptr() + active.n_elem);
}
