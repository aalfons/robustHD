/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#include <R.h>
#include "fastLars.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// ******************************************
// control class definitions for correlations
// ******************************************

// The control classes that handle how the correlations are computed.  They
// store the values for the additional control parameters and have a cor()
// method to call the corresponding correlation function.

// -------------------------------------
// control class for Pearson correlation
// -------------------------------------

class CorPearsonControl {
public:
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// method to compute correlation
double CorPearsonControl::cor(const vec& x, const vec& y) {
	return corPearson(x, y);
}


// ----------------------------------------------------------
// control class for Huber correlation based on winsorization
// ----------------------------------------------------------

// ATTENTION: the data are assumed to be standardized (this is done in R)

class CorHuberControl {
public:
	double c;
	double prob;
	double tol;
	// constructors
	CorHuberControl(double &, double&, double&);
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// constructors
inline CorHuberControl::CorHuberControl(double & _c, double& _p, double& _t) {
	c = _c;
	prob = _p;
	tol = _t;
}

// method to compute correlation
double CorHuberControl::cor(const vec& x, const vec& y) {
	return corHuberBi(x, y, c, prob, tol);
}


// ***************************************************
// (robust) least angle regression variable sequencing
// ***************************************************

// variable sequencing via (robust) least angle regression
// Armadillo library is used for linear algebra
// ATTENTION: the data are assumed to be standardized (this is done in R)
// x ............ predictor matrix
// y ............ response
// sMax ......... number of predictors to be sequenced
// corControl ... control object to compute correlation
// robust ....... check robust correlation matrix for positive definiteness?
// scaleFun ..... R function to compute scale estimates
// ncores ....... number of processor cores for parallel computing
// parallel computing is only used for expensive computations with all
// inactive predictors, otherwise there is no speedup due to overhead
template <class CorControl>
uvec fastLars(const mat& x, const vec& y, const uword& sMax,
		CorControl corControl, const bool& robust, SEXP scaleFun,
		int& ncores) {
	// initializations
	const uword p = x.n_cols;

	// STEP 1: find first ranked predictor
	// compute correlations with response
	vec corY(p);
	#pragma omp parallel for num_threads(ncores) schedule(dynamic)
	for(uword j = 0; j < p; j++) {
		corY(j) = corControl.cor(x.unsafe_col(j), y);
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
		vec xx = x.unsafe_col(active(k-1));
		#pragma omp parallel for num_threads(ncores) schedule(dynamic)
		for(uword j = 0; j < m; j++) {
			vec xj = x.unsafe_col(inactive(j));
			R(inactive(j), k-1) = corControl.cor(xj, xx);
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
        	if(robust) {
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
        	}
            // compute quantities necessary for computing the step size
        	mat invG = solve(G, eye<mat>(k, k));
            a = pow(as_scalar(ones<rowvec>(k) * invG * ones<vec>(k)), -0.5);
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
		r.insert_rows(k, 1, false);	// do not initialize new memory
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
SEXP R_fastLars(SEXP R_x, SEXP R_y, SEXP R_sMax, SEXP R_robust, SEXP R_c,
		SEXP R_prob, SEXP R_tol, SEXP scaleFun, SEXP R_ncores) {
	// data initializations
	NumericMatrix Rcpp_x(R_x);						// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);				// reuse memory
	NumericVector Rcpp_y(R_y);			// response
	vec y(Rcpp_y.begin(), n, false);	// reuse memory
	uword sMax = as<uword>(R_sMax);
	bool robust = as<bool>(R_robust);
	int ncores = as<int>(R_ncores);
	// call native C++ function
	uvec active;
	if(robust) {
		double c = as<double>(R_c);
		double prob = as<double>(R_prob);
		double tol = as<double>(R_tol);
		CorHuberControl corControl(c, prob, tol);
		active = fastLars(x, y, sMax, corControl, robust, scaleFun, ncores) + 1;
	} else {
		CorPearsonControl corControl;
		active = fastLars(x, y, sMax, corControl, false, scaleFun, ncores) + 1;
	}
	// call native C++ function
	return wrap(active.memptr(), active.memptr() + active.n_elem);
}
