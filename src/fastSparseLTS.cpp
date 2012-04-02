/*
 * Author: Andreas Alfons
 *         K.U.Leuven
 */

#include "fastSparseLTS.h"

using namespace Rcpp;
using namespace Eigen;
using namespace std;


// compute sign of a numeric value
int sign(const double& x) {
	return (x > 0) - (x < 0);
}

//// create sequence of integers (starting with 0)
//VectorXi seqLen(const int& n) {
//	VectorXi sequence(n);
//	for(int i = 0; i < n; i++) {
//		sequence(i) = i;
//	}
//	return sequence;
//}

int findMaxActive(const int& n, const int& p, const bool& useIntercept) {
	int maxActive = n - useIntercept;
	if(p < maxActive) {
		maxActive = p;
	}
	return maxActive;
}

// find new active predictors
// absCor ... absolute correlations of inactive variables with current response
// maxCor ... maximum absolute correlation
VectorXi findNewActive(const VectorXd& absCor, const double& maxCor) {
	const int m = absCor.size();	// number of inactive predictors
	int k = 0;						// number of new active predictors
	VectorXi newActive(m);
	for(int j = 0; j < m; j++) {
		if(absCor(j) >= maxCor) {
			newActive(k) = j;
			k++;
		}
	}
	return newActive.head(k);
}

// compute step size in the direction of the equiangular vector
// corActiveY.. ... correlations of active variables with current response
// corInactiveY ... correlations of inactive variables with current response
// corActiveU ..... correlations of active variables with equiangular vector
// corInactiveU ... correlations of inactive variables with equiangular vector
// eps ............ small numerical value (effective zero)
double findStep(const double& corActiveY, const VectorXd& corInactiveY,
		const double& corActiveU, const VectorXd& corInactiveU,
		const double& eps) {
	const int m = corInactiveY.size();	// number of inactive variables
	int k = 0;							// number of positive values
	double step;						// step size
	VectorXd steps(2*m);				// vector of all values to consider
	// first direction
	for(int j = 0; j < m; j++) {
		step = (corActiveY - corInactiveY(j))/(corActiveU - corInactiveU(j));
		if(step > eps) {
			steps(k) = step;
			k++;
		}
	}
	// second direction
	for(int j = 0; j < m; j++) {
		step = (corActiveY + corInactiveY(j))/(corActiveU + corInactiveU(j));
		if(step > eps) {
			steps(k) = step;
			k++;
		}
	}
	// find and return step size
	step = corActiveY/corActiveU;	// maximum possible step;
	if(k > 0) {
		double smallestPositive = steps.head(k).minCoeff();	// smallest positive value
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
VectorXi findDrops(const VectorXd& beta, const VectorXi& active,
		const VectorXd& w, const double& eps, double& step) {
	const int k = active.size();	// number of active variables
	int s = 0;	// number of active variables with sign change
	int d = 0;	// number of variables to be dropped
	VectorXd steps(k);
	VectorXi drops(k);
	double tmp;
	// for each variable, compute step size where sign change would take place,
	// and keep track of indices of variables that are potentially dropped
	for(int j = 0; j < k; j++) {
		tmp = -beta(active(j))/w(j);
		if(tmp > eps) {
			steps(s) = tmp;
			drops(s) = j;
			s++;
		}
	}
	if(s > 0) {
		// check if sign change occurs before the designated step size
		// if so, adjust step size and find variables to be dropped
		double smallestPositive = steps.head(s).minCoeff();	// smallest positive value
		if(smallestPositive < step) {
			step = smallestPositive;
			for(int j = 0; j < s; j++) {
				if(steps(j) == smallestPositive) {
					drops(d) = drops(j);
					d++;
				}
			}
		}
	}
	// if there are no sign changes or sign change would occur after the
	// designated step size, an empty vector is returned
	return drops.head(d);
}

// define new data type for sorting
// .first component is index
// .second component is value
typedef pair<int, double> sortData;

bool sortDataLess(const sortData& left, const sortData& right) {
	return left.second < right.second;
}

// find indices of h smallest observations
// those are not ordered, but they are the h smallest
VectorXi findSmallest(const VectorXd& x, const int& h) {
	// initialize data structure for sorting
	const int n = x.size();
	vector<sortData> foo(n);
	for(int i = 0; i < n; i++) {
		foo[i] = sortData(i, x(i));
	}
	// call STL's nth_element()
	nth_element(foo.begin(), foo.begin()+h, foo.end(), sortDataLess);
	// construct and return vector of indices
	VectorXi indices(h);
	for(int i = 0; i < h; i++) {
		sortData bar = foo[i];
		indices(i) = bar.first;
	}
	return indices;
}

// R interface to findSmallest()
SEXP R_findSmallest(SEXP R_x, SEXP R_h) {
	NumericVector Rcpp_x(R_x);
	Map<VectorXd> x(Rcpp_x.begin(), Rcpp_x.size());	// reuse memory
	int h = as<int>(R_h);
	VectorXi indices = findSmallest(x, h);	// call native C++ function
	return wrap(indices);
}

//// find indices of h smallest observations
//VectorXi partialOrder(const VectorXd& x, const int& h) {
//	// initialize data structure for sorting
//	const int n = x.size();
//	vector<sortData> foo(n);
//	for(int i = 0; i < n; i++) {
//		foo[i] = sortData(i, x(i));
//	}
//	// call STL's partial_sort()
//	partial_sort(foo.begin(), foo.begin()+h, foo.end(), sortDataLess);
//	// construct and return vector of indices
//	VectorXi indices(h);
//	for(int i = 0; i < h; i++) {
//		sortData bar = foo[i];
//		indices(i) = bar.first;
//	}
//	return indices;
//}
//
//// R interface to partialOrder()
//SEXP R_partialOrder(SEXP R_x, SEXP R_h) {
//	NumericVector Rcpp_x(R_x);
//	Map<VectorXd> x(Rcpp_x.begin(), Rcpp_x.size());	// reuse memory
//	int h = as<int>(R_h);
//	VectorXi indices = partialOrder(x, h);	// call native C++ function
//	return wrap(indices);
//}
//
//// R interface to partial sort of h smallest observations
//SEXP R_partialSort(SEXP R_x, SEXP R_h) {
//	NumericVector Rcpp_x(R_x);
//	NumericVector x(Rcpp_x.begin(), Rcpp_x.end());	// copy original vector
//	int h = as<int>(R_h);
//	// call STL's partial_sort()
//	partial_sort(x.begin(), x.begin()+h, x.end());
//	return x;
//}

// order observations
VectorXi order(const VectorXd& x) {
	// initialize data structure for sorting
	const int n = x.size();
	vector<sortData> foo(n);
	for(int i = 0; i < n; i++) {
		foo[i] = sortData(i, x(i));
	}
	// call STL's sort()
	sort(foo.begin(), foo.end(), sortDataLess);
	// construct and return vector of indices
	VectorXi indices(n);
	for(int i = 0; i < n; i++) {
		sortData bar = foo[i];
		indices(i) = bar.first;
	}
	return indices;
}

VectorXi bestSubsets(const VectorXd& x, const MatrixXi& subsets, int& nkeep) {
	// initializations
	int n = x.size(), h = subsets.rows();
	VectorXi ord = order(x), keep(nkeep);
	keep(0) = ord(0);	// subset with lowest value of objective function
	// loop over subsets until nkeep unique subsets are found
	int j = 1, k = 1;
	while((k < nkeep) && (j < n)) {
		if(x(keep(k-1)) < x(ord(j))) {
			// keep the current subset if its value for objective function is
			// larger than for the previous subset
			keep(k) = ord(j);
			k++;
		} else {
			// equal value for objective function as for the previous subset
			// sort indices of previous and current subsets and check if they
			// are equal
			VectorXi previous = subsets.col(keep(k-1));
			sort(&previous.coeffRef(0), &previous.coeffRef(0)+h);
			VectorXi current = subsets.col(ord(j));
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
				keep(k) = ord(j);
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


// barebones version of the lasso for fixed lambda
// Eigen library is used for linear algebra
// ATTENTION: intercept is returned through corresponding parameter
// x .............. predictor matrix
// y .............. response
// lambda ......... penalty parameter
// useSubset ...... logical indicating whether lasso should be computed on a
//                  subset
// subset ......... indices of subset on which lasso should be computed
// useIntercept ... logical indicating whether intercept should be included
// intercept ...... intercept is returned through this parameter
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
VectorXd fastLasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useSubset, const VectorXi& subset, const bool& useIntercept,
		double& intercept, const double& eps, const bool& useGram) {

	// data initializations
	int n, p = x.cols();
	MatrixXd xs;
	VectorXd ys;
	if(useSubset) {
		n = subset.size();
		xs.resize(n, p);
		ys.resize(n);
		int s;
		for(int i = 0; i < n; i++) {
			s = subset(i);
			xs.row(i) = x.row(s);
			ys(i) = y(s);
		}
	} else {
		n = x.rows();
		xs = x;	// does this copy memory?
		ys = y;	// does this copy memory?
	}
	double rescaledLambda = n * lambda / 2;

	// center data and store means
	RowVectorXd meanX;
	double meanY;
	if(useIntercept) {
		meanX = xs.colwise().mean();	// columnwise means of predictors
		xs.rowwise() -= meanX;			// sweep out columnwise means
		meanY = ys.mean();				// mean of response
		for(int i = 0; i < n; i++) {
			ys(i) -= meanY;				// sweep out mean
		}
	} else {
		meanY = 0;		// just to avoid warning, this is never used
//		intercept = 0;	// zero intercept
	}

	// some initializations
	VectorXi inactive(p);	// inactive predictors
	int m = 0;				// number of inactive predictors
	VectorXi ignores;		// indicates variables to be ignored
	int s = 0;				// number of ignored variables

	// normalize predictors and store norms
	RowVectorXd normX = xs.colwise().norm();	// columnwise norms
	double epsNorm = eps * sqrt(n);	// R package 'lars' uses n, not n-1
	for(int j = 0; j < p; j++) {
		if(normX(j) < epsNorm) {
			// variance is too small: ignore variable
			ignores.append(j, s);
			s++;
			// set norm to tolerance to avoid numerical problems
			normX(j) = epsNorm;
		} else {
			inactive(m) = j;		// add variable to inactive set
			m++;					// increase number of inactive variables
		}
		xs.col(j) /= normX(j);		// sweep out norm
	}
	// resize inactive set and update number of variables if necessary
	if(m < p) {
		inactive.conservativeResize(m);
		p = m;
	}

	// compute Gram matrix if requested (saves time if number of variables is
	// not too large)
	MatrixXd Gram;
	if(useGram) {
		Gram.noalias() = xs.transpose() * xs;
	}

	// further initializations for iterative steps
	RowVectorXd corY;
	corY.noalias() = ys.transpose() * xs;	// current correlations
	VectorXi active;	// active predictors
	int k = 0;			// number of active predictors
	VectorXd previousBeta = VectorXd::Zero(p+s), currentBeta = VectorXd::Zero(p+s);	// previous and current regression coefficients
	double previousLambda = 0, currentLambda = 0;	// previous and current penalty parameter
	VectorXi drops;		// indicates variables to be dropped
	VectorXd signs;		// keep track of sign of correlations for the active variables (double precision is necessary for solve())
	MatrixXd L;			// Cholesky L of Gram matrix of active variables
	int rank = 0;		// rank of Cholesky L
	int maxActive = findMaxActive(n, p, useIntercept);	// maximum number of variables to be sequenced

	// modified LARS algorithm for lasso solution
	while(k < maxActive) {

		// extract current correlations of inactive variables
		VectorXd corInactiveY(m);
		for(int j = 0; j < m; j++) {
			corInactiveY(j) = corY(inactive(j));
		}
		// compute absolute values of correlations and find maximum
		VectorXd absCorInactiveY = corInactiveY.cwiseAbs();
		double maxCor = absCorInactiveY.maxCoeff();
		// update current lambda
		if(k == 0) {	// no active variables
			previousLambda = maxCor;
		} else {
			previousLambda = currentLambda;
		}
		currentLambda = maxCor;
		if(currentLambda <= rescaledLambda) break;

		if(drops.size() == 0) {
			// new active variables
			VectorXi newActive = findNewActive(absCorInactiveY, maxCor - eps);
			// do calculations for new active variables
			for(int j = 0; j < newActive.size(); j++) {
				// update Cholesky L of Gram matrix of active variables
				// TODO: put this into void function
				int newJ = inactive(newActive(j));
				VectorXd xNewJ;
				double newX;
				if(useGram) {
					newX = Gram(newJ, newJ);
				} else {
					xNewJ = xs.col(newJ);
					newX = xNewJ.squaredNorm();
				}
				double normNewX = sqrt(newX);
				if(k == 0) {	// no active variables, L is empty
					L.resize(1,1);
					L(0, 0) = normNewX;
					rank = 1;
				} else {
					VectorXd oldX(k);
					if(useGram) {
						for(int j = 0; j < k; j++) {
							oldX(j) = Gram(active(j), newJ);
						}
					} else {
						for(int j = 0; j < k; j++) {
							oldX(j) = xNewJ.dot(xs.col(active(j)));
						}
					}
					VectorXd l = L.triangularView<Lower>().solve(oldX);
					double lkk = newX - l.squaredNorm();
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
						L.conservativeResize(k+1, k+1);
						for(int j = 0; j < k; j++) {
							L(k, j) = l(j);
							L(j, k) = l(j);
						}
						L(k,k) = lkk;
					}
				}
				// add new variable to active set or drop it for good
				// in case of singularity
				if(rank == k) {
					// singularity: drop variable for good
					ignores.append(newJ, s);
					s++;	// increase number of ignored variables
					p--;	// decrease number of variables
					if(p < maxActive) {
						// adjust maximum number of active variables
						maxActive = p;
					}
				} else {
					// no singularity: add variable to active set
					active.append(newJ, k);
					// keep track of sign of correlation for new active variable
					signs.append(sign(corY(newJ)), k);
					k++;	// increase number of active variables
				}
			}
			// remove new active or ignored variables from inactive variables
			// and corresponding vector of current correlations
			inactive.remove(newActive);
			corInactiveY.remove(newActive);
			m = inactive.size();	// update number of inactive variables
		}
		// prepare for computation of step size
		// here double precision of signs is necessary
		VectorXd b = L.triangularView<Lower>().solve(signs);
		VectorXd G = L.triangularView<Upper>().solve(b);
		// correlations of active variables with equiangular vector
		double corActiveU = 1/sqrt(G.dot(signs));
		// coefficients of active variables in linear combination forming the
		// equiangular vector
		VectorXd w = G * corActiveU;	// note that this has the right signs
		// equiangular vector
		VectorXd u;
		if(!useGram) {
			// we only need equiangular vector if we don't use the precomputed
			// Gram matrix, otherwise we can compute the correlations directly
			// from the Gram matrix
			u = VectorXd::Zero(n);
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < k; j++) {
					u(i) += xs(i, active(j)) * w(j);
				}
			}
		}
		// compute step size in equiangular direction
		double step;
		if(k < maxActive) {
			// correlations of inactive variables with equiangular vector
			VectorXd corInactiveU(m);
			if(useGram) {
				for(int j = 0; j < m; j++) {
					corInactiveU(j) = 0;
					for(int i = 0; i < k; i++) {
						corInactiveU(j) += w(i) * Gram(active(i), inactive(j));
					}
				}
			} else {
				for(int j = 0; j < m; j++) {
					corInactiveU(j) = u.dot(xs.col(inactive(j)));
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
		for(int j = 0; j < k; j++) {
			currentBeta(active(j)) += step * w(j);
		}
		// update current correlations
		if(useGram) {
			// we also need to do this for active variables, since they may be
			// dropped at a later stage
			// TODO: computing a vector step * w in advance may save some computation time
			for(int j = 0; j < Gram.cols(); j++) {
				for(int i = 0; i < k; i++) {
					corY(j) -= step * w(i) * Gram(active(i), j);
				}
			}
		} else {
			ys -= step * u;	// take step in equiangular direction
			corY.noalias() = ys.transpose() * xs;	// update correlations
		}
		// drop variables if necessary
		if(drops.size() > 0) {
			// downdate Cholesky L
			// TODO: put this into void function
			for(int j = drops.size()-1; j >= 0; j--) {
				// variables need to be dropped in descending order
				int drop = drops(j);	// index with respect to active set
				// modify upper triangular part as in R package 'lars'
				// simply deleting columns is not enough, other computations
				// necessary but complicated due to Fortran code
				L.removeCol(drop);
				VectorXd z = VectorXd::Constant(k, 1, 1);
				k--;	// decrease number of active variables
				for(int i = drop; i < k; i++) {
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
						for(int j = i+1; j < k; j++) {
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
				// TODO: removing all rows together may save some computation time
				L.conservativeResize(k, NoChange);
				rank--;
			}
			// mirror lower triangular part
			for(int j = 0; j < k; j++) {
				for(int i = j+1; i < k; i++) {
					L(i,j) = L(j,i);
				}
			}
			// add dropped variables to inactive set and make sure
			// coefficients are 0
			inactive.conservativeResize(m + drops.size());
			for(int j = 0; j < drops.size(); j++) {
				int newInactive = active(drops(j));
				inactive(m + j) = newInactive;
				currentBeta(newInactive) = 0;	// make sure coefficient is 0
			}
			m = inactive.size();	// update number of inactive variables
			// drop variables from active set and sign vector
			// number of active variables is already updated above
			active.remove(drops);
			signs.remove(drops);
		}
	}

	// interpolate coefficients for given lambda
	p = currentBeta.size();	// reset number of variables to include ignored ones
	VectorXd beta(p);		// final coefficient vector
	if(rescaledLambda <= currentLambda) {
		beta = currentBeta;
	} else if(rescaledLambda >= previousLambda) {
		// strict inequality only possible if lambda larger than largest lambda
		beta = previousBeta;	// includes case of equality
	} else {
		// penalty parameter strictly within two steps: interpolate coefficients
		beta = ((rescaledLambda - currentLambda) * previousBeta +
				(previousLambda - rescaledLambda) * currentBeta) /
				(previousLambda - currentLambda);
	}

	// transform coefficients back
	for(int j = 0; j < p; j++) {
		beta(j) /= normX(j);
	}
	if(useIntercept) {
		intercept = meanY - beta.dot(meanX);
	}

	return beta;
}

// wrapper function used internally by fastSparseLTS(), which returns residuals
// and value of objective function through corresponding parameters
VectorXd fastLasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const VectorXi& subset, const bool& useIntercept, double& intercept,
		VectorXd& residuals, double& crit, const double& eps,
		const bool& useGram) {
	// compute coefficients
	VectorXd coefficients = fastLasso(x, y, lambda, true, subset,
			useIntercept, intercept, eps, useGram);
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

// wrapper function used for R interface, which returns fitted values and
// residuals through corresponding parameters
VectorXd fastLasso(const MatrixXd& x, const VectorXd& y, const double& lambda,
		const bool& useSubset, const VectorXi& subset, const bool& useIntercept,
		double& intercept, VectorXd& fitted, VectorXd& residuals,
		const double& eps, const bool& useGram) {
	// compute coefficients
	VectorXd coefficients = fastLasso(x, y, lambda, useSubset, subset,
			useIntercept, intercept, eps, useGram);
	// compute fitted values
	fitted.noalias() = x * coefficients;
	if(useIntercept) {
		for(int i = 0; i < fitted.size(); i++) {
			fitted(i) += intercept;
		}
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
	NumericMatrix Rcpp_x(R_x);	// predictor matrix
	const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
	Map<MatrixXd> x(Rcpp_x.begin(), n, p);	// reuse memory
	NumericVector Rcpp_y(R_y);	// response
	Map<VectorXd> y(Rcpp_y.begin(), n);		// reuse memory
	double lambda = as<double>(R_lambda);
	bool useSubset = as<bool>(R_useSubset);
	IntegerVector Rcpp_subset(R_subset);	// subset to use for computation
	Map<VectorXi> subset(Rcpp_subset.begin(), Rcpp_subset.size());
	bool useIntercept = as<bool>(R_intercept);
	double intercept;
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function and return results as list
	VectorXd fitted, residuals;
	VectorXd coefficients = fastLasso(x, y, lambda, useSubset, subset,
			useIntercept, intercept, fitted, residuals, eps, useGram);
	NumericVector R_coefficients = wrap(coefficients);
	if(useIntercept) {
		R_coefficients.push_front(intercept);	// prepend intercept
	}
	return List::create(
			Named("coefficients") = R_coefficients,
			Named("fitted.values") = fitted,
			Named("residuals") = residuals
			);
}


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
				useIntercept, intercept, eps, useGram);
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


// C-step for sparse LTS for k-th subset
// Eigen library is used for linear algebra
// ATTENTION: coefficients, intercept, residuals and subset are returned
//            through corresponding parameters
// x .............. predictor matrix
// y .............. response
// lambda ......... penalty parameter
// k .............. index of subset to consider
// residuals ...... matrix of residuals for all subsets
// h .............. subset size
// subsets ........ matrix of all subsets
// coefficients ... matrix of regression coefficients for all subsets
// useIntercept ... logical indicating whether intercept should be included
// intercepts ..... vector of intercepts for all subsets (if included in the model)
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
VectorXd fastCStep(const MatrixXd& x, const VectorXd& y, const double& lambda,
		VectorXd& residuals, const int& h, VectorXi& subset,
		const bool& useIntercept, double& intercept, double& crit,
		const double& eps, const bool& useGram) {
	// update subset
	subset = findSmallest(residuals.cwiseAbs(), h);
	// compute lasso solution for new subset
	return fastLasso(x, y, lambda, subset, useIntercept,
			intercept, residuals, crit, eps, useGram);
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
	double intercept;
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function and return results as list
	VectorXi subset;
	double crit;
	VectorXd coefficients = fastCStep(x, y, lambda, residuals, h, subset,
			useIntercept, intercept, crit, eps, useGram);
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


// sparse least trimmed squares
// Eigen library is used for linear algebra
// x .............. predictor matrix
// y .............. response
// lambda ......... penalty parameter
// subsets ........ matrix of initial subsets
// useIntercept ... logical indicating whether intercept should be included
// intercept ...... intercept (if included in the model)
// tol ............ numerical tolerance for convergence
// eps ............ small numerical value (effective zero)
// useGram ........ logical indicating whether Gram matrix should be computed
//                  in advance
VectorXd fastSparseLTS(const MatrixXd& x, const VectorXd& y,
		const double& lambda, Map<MatrixXi>& subsets, const bool& useIntercept,
		double& intercept, int& nkeep, const int& ncstep, VectorXi& best,
		VectorXd& residuals, double& crit, double& center, double& scale,
		const double& tol, const double& eps, const bool& useGram) {
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
				intercept, residuals, crit, eps, useGram);
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
					useIntercept, intercept, crit, eps, useGram);
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
					useIntercept, intercept, crit, eps, useGram);
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
	best = subsets.col(which);
	residuals = residMat.col(which);
	crit = crits(which);
	center = subsetMean(residuals, best);  		// residual center
	scale = partialScale(residuals, center, h);	// uncorrected residual scale
	if(useIntercept) {
		intercept = intercepts(which);
	}
	return coefMat.col(which);
}

// R interface to fastSparseLTS()
// initial subsets are constructed in R and passed down to C++
SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_subsets,
		SEXP R_intercept, SEXP R_nkeep, SEXP R_ncstep, SEXP R_tol, SEXP R_eps,
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
	double intercept;
	int nkeep = as<int>(R_nkeep);
	int ncstep = as<int>(R_ncstep);
	double tol = as<double>(R_tol);
	double eps = as<double>(R_eps);
	bool useGram = as<bool>(R_useGram);
	// call native C++ function
	VectorXi best;
	VectorXd residuals;
	double crit, center, scale;
	VectorXd coefficients = fastSparseLTS(x, y, lambda, subsets, useIntercept,
			intercept, nkeep, ncstep, best, residuals, crit, center, scale,
			tol, eps, useGram);
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
