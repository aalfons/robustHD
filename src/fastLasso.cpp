/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#include "fastLasso.h"

using namespace Rcpp;
using namespace arma;


// find maximum number of active variables
// n .............. number of observations
// p .............. number of predictors
// useIntercept ... logical indicating whether model has an intercept
uword findMaxActive(const uword& n, const uword& p, const bool& useIntercept) {
	uword maxActive = n - useIntercept;
	if(p < maxActive) {
		maxActive = p;
	}
	return maxActive;
}

// compute step size in the direction of the equiangular vector
// corActiveY.. ... correlations of active variables with current response
// corInactiveY ... correlations of inactive variables with current response
// corActiveU ..... correlations of active variables with equiangular vector
// corInactiveU ... correlations of inactive variables with equiangular vector
// eps ............ small numerical value (effective zero)
double findStep(const double& corActiveY, const vec& corInactiveY,
		const double& corActiveU, const vec& corInactiveU,
		const double& eps) {
	// construct vector of all values to consider
	vec steps = join_cols((corActiveY - corInactiveY)/(corActiveU - corInactiveU),
			(corActiveY + corInactiveY)/(corActiveU + corInactiveU));
	steps = steps.elem(find(steps > eps));
	// find and return step size
	double step = corActiveY/corActiveU;      // maximum possible step;
	if(steps.n_elem > 0) {
		double smallestPositive = steps.min();  // smallest positive value
		if(smallestPositive < step) {
			step = smallestPositive;
		}
	}
	return step;
}

// adjust step size if any sign changes before the designated step size,
// and return the corresponding variables to be dropped
// beta   ... current regression coefficients
// active ... indices of inactive variables
// w ........ coefficients of active variables in linear combination forming
//            the equiangular vector
// eps ...... small numerical value (effective zero)
// step ..... step size in direction of equiangular vector
uvec findDrops(const vec& beta, const uvec& active, const vec& w,
		const double& eps, double& step) {
	// for each variable, compute step size where sign change would take place,
	// and keep track of indices of variables that are potentially dropped
	vec steps = -beta.elem(active) / w;
	uvec drops = find(steps > eps);
	if(drops.n_elem > 0) {
		// check if sign change occurs before the designated step size
		// if so, adjust step size and find variables to be dropped
		steps = steps.elem(drops);
		double smallestPositive = steps.min();
		if(smallestPositive < step) {
			step = smallestPositive;
			drops = drops.elem(find(steps == smallestPositive));
		} else drops.reset();
	}
	// if there are no sign changes or sign change would occur after the
	// designated step size, an empty vector is returned
	return drops;
}


// barebones version of the lasso for fixed lambda
// Armadillo library is used for linear algebra
// ATTENTION: intercept is returned through corresponding parameter
// x .............. predictor matrix
// y .............. response
// lambda ......... penalty parameter
// useSubset ...... logical indicating whether lasso should be computed on a
//                  subset
// subset ......... indices of subset on which lasso should be computed
// useIntercept ... logical indicating whether intercept should be included
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
// intercept ...... intercept is returned through this parameter
vec fastLasso(const mat& x, const vec& y, const double& lambda,
		const bool& useSubset, const uvec& subset, const bool& useIntercept,
		const double& eps, const bool& useGram, double& intercept) {

	// data initializations
	uword n, p = x.n_cols;
	mat xs;
	vec ys;
	if(useSubset) {
		n = subset.n_elem;
		xs.set_size(n, p);
		ys.set_size(n);
		uword s;
		for(uword i = 0; i < n; i++) {
			s = subset(i);
			xs.row(i) = x.row(s);
			ys(i) = y(s);
		}
	} else {
		n = x.n_rows;
		xs = x;	// does this copy memory?
		ys = y;	// does this copy memory?
	}
	double rescaledLambda = n * lambda / 2;

	// center data and store means
	rowvec meanX;
	double meanY;
	if(useIntercept) {
		meanX = mean(xs, 0);		  // columnwise means of predictors
		for(uword j = 0; j < p; j++) {
			xs.col(j) -= meanX(j);	// sweep out columnwise means
		}
		meanY = mean(ys);	// mean of response
		ys -= meanY;		  // sweep out mean
	} else {
		meanY = 0;		  // just to avoid warning, this is never used
//		intercept = 0;	// zero intercept
	}

	// compute norms and find variables with too small a norm
	uvec inactive = seqLen(p);
	rowvec normX = sqrt(sum(xs % xs, 0));	// columnwise norms
	double epsNorm = eps * sqrt(n);	      // R package 'lars' uses n, not n-1
	uvec ignores = find(normX < epsNorm);	// indicates variables to be ignored
	uword s = ignores.n_elem;
	for(sword j = s-1; j >= 0; j--) { // reverse order (requires signed integer)
		uword i = ignores(j);
		// set norm to tolerance to avoid numerical problems
		normX(i) = epsNorm;
		// remove ignored variable from inactive set (hence reverse order)
		inactive.shed_row(i);
	}
	// normalize predictors
	for(uword j = 0; j < p; j++) {
		xs.col(j) /= normX(j);		// sweep out norm
	}

	// compute Gram matrix if requested (saves time if number of variables is
	// not too large)
	mat Gram;
	if(useGram) {
		Gram = trans(xs) * xs;
	}

	// further initializations
	rowvec corY = conv_to<rowvec>::from(ys) * xs;	// current correlations
  vec beta(p);                // final coefficients
  uword m = inactive.n_elem;  // number of inactive predictors
  if(m < p) p = m;            // update number of predictors if necessary

  // compute lasso solution
  if(p == 1) {

    // special case of only one variable (with sufficiently large norm)
    uword j = inactive(0);          
    // set maximum step size in the direction of that variable
    double maxStep = corY(j);
    if(maxStep < 0) maxStep = -maxStep; // absolute value
    // compute coefficients for least squares solution
    vec betaLS = solve(xs.unsafe_col(j), ys);
    // compute lasso coefficients
    beta.zeros();
    if(rescaledLambda < maxStep) {
      // interpolate towards least squares solution
      beta(j) = (maxStep - rescaledLambda) * betaLS(0) / maxStep;
    }

  } else {

    // further initializations for iterative steps
    uvec active;  // active predictors
    uword k = 0;	// number of active predictors
    // previous and current regression coefficients
  	vec previousBeta = zeros(p+s), currentBeta = zeros(p+s);
  	// previous and current penalty parameter
    double previousLambda = R_PosInf, currentLambda = R_PosInf;
  	// indicates variables to be dropped
    uvec drops;
  	// keep track of sign of correlations for the active variables 
    // (double precision is necessary for solve())
    vec signs;
  	// Cholesky L of Gram matrix of active variables
    mat L;
  	uword rank = 0;		// rank of Cholesky L
    // maximum number of variables to be sequenced
  	uword maxActive = findMaxActive(n, p, useIntercept);

    // modified LARS algorithm for lasso solution
    while((k < maxActive)) {

  		// extract current correlations of inactive variables
  		vec corInactiveY = corY.elem(inactive);
  		// compute absolute values of correlations and find maximum
  		vec absCorInactiveY = abs(corInactiveY);
  		double maxCor = absCorInactiveY.max();
  		// update current lambda
  		if(k == 0) {	// no active variables
  			previousLambda = maxCor;
  		} else {
  			previousLambda = currentLambda;
  		}
  		currentLambda = maxCor;
  		if(currentLambda <= rescaledLambda) break;

  		if(drops.n_elem == 0) {
  			// new active variables
  			uvec newActive = find(absCorInactiveY >= (maxCor - eps));
  			// do calculations for new active variables
  			for(uword j = 0; j < newActive.n_elem; j++) {
  				// update Cholesky L of Gram matrix of active variables
  				// this cannot be put into its own void function since
  				// insert_rows() doesn't work with referenced matrices
  				uword newJ = inactive(newActive(j));
  				vec xNewJ;
  				double newX;
  				if(useGram) {
  					xNewJ = Gram.unsafe_col(newJ);	// reuses memory
  					newX = xNewJ(newJ);
  				} else {
  					xNewJ = xs.unsafe_col(newJ);	  // reuses memory
  					newX = accu(xNewJ % xNewJ);
  				}
  				double normNewX = sqrt(newX);
  				if(k == 0) {	// no active variables, L is empty
  					L.set_size(1,1);
  					L(0, 0) = normNewX;
  					rank = 1;
  				} else {
  					vec oldX;
  					if(useGram) {
  						oldX = xNewJ.elem(active);
  					} else {
  						oldX.set_size(k);
  						for(uword j = 0; j < k; j++) {
  							oldX(j) = dot(xNewJ, xs.unsafe_col(active(j)));
  						}
  					}
  					vec l = solve(trimatl(L), oldX);
  					double lkk = newX - accu(l % l);
  					// check if L is machine singular
  					if(lkk > eps) {
  						// no singularity: update Cholesky L
  						lkk = sqrt(lkk);
  						rank++;
  						// add new row and column to Cholesky L
  						// this is a little trick: sometimes we need
  						// lower triangular matrix, sometimes upper
  						// hence we define quadratic matrix and use
  						// triangularView() to interpret matrix the
  						// correct way
  						// insert row and column without initializing memory
  						// (set_size() and reshape() have strange behavior)
  						L.insert_rows(k, 1, false);
  						L.insert_cols(k, 1, false);
  						// fill new parts of the matrix
  						for(uword j = 0; j < k; j++) {
  							L(k, j) = l(j);
  							L(j, k) = l(j);
  						}
  						L(k, k) = lkk;
  					}
  				}
  				// add new variable to active set or drop it for good
  				// in case of singularity
  				if(rank == k) {
  					// singularity: drop variable for good
  					ignores.insert_rows(s, 1, false);	// do not initialize new memory
  					ignores(s) = newJ;
  					s++;	// increase number of ignored variables
  					p--;	// decrease number of variables
  					if(p < maxActive) {
  						// adjust maximum number of active variables
  						maxActive = p;
  					}
  				} else {
  					// no singularity: add variable to active set
  					active.insert_rows(k, 1, false);	// do not initialize new memory
  					active(k) = newJ;
  					// keep track of sign of correlation for new active variable
  					signs.insert_rows(k, 1, false);		// do not initialize new memory
  					signs(k) = sign(corY(newJ));
  					k++;	// increase number of active variables
  				}
  			}
  			// remove new active or ignored variables from inactive variables
  			// and corresponding vector of current correlations
  			for(sword j = newActive.n_elem-1; j >= 0; j--) {	// reverse order
  				uword i = newActive(j);
  				inactive.shed_row(i);
  				corInactiveY.shed_row(i);
  			}
  			m = inactive.n_elem;	// update number of inactive variables
  		}
  		// prepare for computation of step size
  		// here double precision of signs is necessary
  		vec b = solve(trimatl(L), signs);
  		vec G = solve(trimatu(L), b);
  		// correlations of active variables with equiangular vector
  		double corActiveU = 1/sqrt(dot(G, signs));
  		// coefficients of active variables in linear combination forming the
  		// equiangular vector
  		vec w = G * corActiveU;	// note that this has the right signs
  		// equiangular vector
  		vec u;
  		if(!useGram) {
  			// we only need equiangular vector if we don't use the precomputed
  			// Gram matrix, otherwise we can compute the correlations directly
  			// from the Gram matrix
  			u = zeros<vec>(n);
  			for(uword i = 0; i < n; i++) {
  				for(uword j = 0; j < k; j++) {
  					u(i) += xs(i, active(j)) * w(j);
  				}
  			}
  		}
  		// compute step size in equiangular direction
  		double step;
  		if(k < maxActive) {
  			// correlations of inactive variables with equiangular vector
  			vec corInactiveU(m);
  			if(useGram) {
  				for(uword j = 0; j < m; j++) {
  					vec gram = Gram.unsafe_col(inactive(j));
  					corInactiveU(j) = dot(w, gram.elem(active));
  				}
  			} else {
  				for(uword j = 0; j < m; j++) {
  					corInactiveU(j) = dot(u, xs.unsafe_col(inactive(j)));
  				}
  			}
  			// compute step size in the direction of the equiangular vector
  			step = findStep(maxCor, corInactiveY, corActiveU, corInactiveU, eps);
  		} else {
  			// last step: take maximum possible step
  			step = maxCor/corActiveU;
  		}
  		// adjust step size if any sign changes and drop corresponding variables
  		drops = findDrops(currentBeta, active, w, eps, step);
  		// update current regression coefficients
  		previousBeta = currentBeta;
  		currentBeta.elem(active) += step * w;
  		// update current correlations
  		if(useGram) {
  			// we also need to do this for active variables, since they may be
  			// dropped at a later stage
  			for(uword j = 0; j < Gram.n_cols; j++) {
  				vec gram = Gram.unsafe_col(j);
  				corY(j) -= step * dot(w, gram.elem(active));
  			}
  		} else {
  			ys -= step * u;	// take step in equiangular direction
  			corY = conv_to<rowvec>::from(ys) * xs;	// might be faster than trans()
  		}
  		// drop variables if necessary
  		if(drops.n_elem > 0) {
  			// downdate Cholesky L
  			// this cannot be put into its own void function since
  			// shed_col() and shed_row() don't work with referenced matrices
  			for(sword j = drops.n_elem-1; j >= 0; j--) {	// reverse order
  				// variables need to be dropped in descending order
  				uword drop = drops(j);	// index with respect to active set
  				// modify upper triangular part as in R package 'lars'
  				// simply deleting columns is not enough, other computations
  				// necessary but complicated due to Fortran code
  				L.shed_col(drop);
  				vec z = ones<vec>(k);
  				k--;	// decrease number of active variables
  				for(uword i = drop; i < k; i++) {
  					double a = L(i,i), b = L(i+1,i);
  					if(b != 0.0) {
  						// compute the rotation
  						double tau, s, c;
  						if(abs(b) > abs(a)) {
  							tau = -a/b;
  							s = 1.0/sqrt(1.0+tau*tau);
  							c = s * tau;
  						} else {
  							tau = -b/a;
  							c = 1.0/sqrt(1.0+tau*tau);
  							s = c * tau;
  						}
  						// update 'L' and 'z';
  						L(i,i) = c*a - s*b;
  						for(uword j = i+1; j < k; j++) {
  							a = L(i,j);
  							b = L(i+1,j);
  							L(i,j) = c*a - s*b;
  							L(i+1,j) = s*a + c*b;
  						}
  						a = z(i);
  						b = z(i+1);
  						z(i) = c*a - s*b;
  						z(i+1) = s*a + c*b;
  					}
  				}
  				L.shed_row(k);
  				rank--;
  			}
  			// mirror lower triangular part
  			L = symmatu(L);
  			// add dropped variables to inactive set and make sure
  			// coefficients are 0
  			inactive.insert_rows(m, drops.n_elem, false);
  			for(uword j = 0; j < drops.n_elem; j++) {
  				uword newInactive = active(drops(j));
  				inactive(m + j) = newInactive;
  				currentBeta(newInactive) = 0;	// make sure coefficient is 0
  			}
  			m = inactive.n_elem;	// update number of inactive variables
  			// drop variables from active set and sign vector
  			// number of active variables is already updated above
  			for(sword j = drops.n_elem-1; j >= 0; j--) {	// reverse order
  				// variables need to be dropped in descending order
  				uword drop = drops(j);	// index with respect to active set
  				// drop variables from active set and sign vector
  				// number of active variables is already updated above
  				active.shed_row(drop);
  				signs.shed_row(drop);
  			}
  		}
  	}

  	// interpolate coefficients for given lambda
    if(k == 0) {
    	// lambda larger than largest lambda from steps of LARS algorithm
  		beta.zeros();
    } else {
    	// penalty parameter within two steps
      if(k == maxActive) {
          // current coefficients are the least squares solution (in the 
          // high-dimensional case, as far along the solution path as possible)
          // current and previous values of the penalty parameter need to be 
          // reset for interpolation
          previousLambda = currentLambda;
          currentLambda = 0;
      }
      // interpolate coefficients
    	beta = ((rescaledLambda - currentLambda) * previousBeta +
  				(previousLambda - rescaledLambda) * currentBeta) /
  				(previousLambda - currentLambda);
    }
  }

	// transform coefficients back
	beta = beta / conv_to<colvec>::from(normX);
	if(useIntercept) {
		intercept = meanY - dot(beta, meanX);
	}

	return beta;
}

// wrapper function used for R interface, which returns fitted values and
// residuals through corresponding parameters
vec fastLasso(const mat& x, const vec& y, const double& lambda,
		const bool& useSubset, const uvec& subset, const bool& useIntercept,
		const double& eps, const bool& useGram, double& intercept, vec& fitted,
		vec& residuals) {
	// compute coefficients
	vec coefficients = fastLasso(x, y, lambda, useSubset, subset,
			useIntercept, eps, useGram, intercept);
	// compute fitted values
	fitted = x * coefficients;
	if(useIntercept) {
		fitted += intercept;
	}
	// compute residuals
	residuals = y - fitted;
	// return coefficients
	return coefficients;
}

// R interface to fastLasso()
SEXP R_fastLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_useSubset,
		SEXP R_subset, SEXP R_intercept, SEXP R_eps, SEXP R_useGram) {
    // data initializations
	NumericMatrix Rcpp_x(R_x);						  // predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	mat x(Rcpp_x.begin(), n, p, false);		  // reuse memory
	NumericVector Rcpp_y(R_y);			        // response
	vec y(Rcpp_y.begin(), n, false);	      // reuse memory
	double lambda = as<double>(R_lambda);
	bool useSubset = as<bool>(R_useSubset);
	uvec subset;
	if(useSubset) {
		IntegerVector Rcpp_subset(R_subset);	// subset to use for computation
		const int h = Rcpp_subset.size();
		subset = uvec(h);
		for(int i = 0; i < h; i++) {
      // can't use the same memory-saving conversion for integer vectors
			subset(i) = Rcpp_subset[i] - 1;
		}
	}
	bool useIntercept = as<bool>(R_intercept);
	double intercept;
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function and return results as list
	vec fitted, residuals;
	vec coefficients = fastLasso(x, y, lambda, useSubset, subset, useIntercept,
			eps, useGram, intercept, fitted, residuals);
	if(useIntercept) {
		// prepend intercept
		coefficients.insert_rows(0, 1, false);
		coefficients(0) = intercept;
	}
	return List::create(
			Named("coefficients") = coefficients,
			Named("fitted.values") = fitted,
			Named("residuals") = residuals
			);
}
