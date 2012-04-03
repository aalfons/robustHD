/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include <R.h>
#include "corHuber.h"

using namespace Rcpp;
using namespace Eigen;


// -------------------
// Pearson correlation
// -------------------

double corPearson(const VectorXd& x, const VectorXd& y) {
	const int n = x.size();
	// compute means
	double meanX = 0, meanY = 0;
	for(int i = 0; i < n; i++) {
		meanX += x(i);
		meanY += y(i);
	}
	meanX /= n;
	meanY /= n;
	// compute Pearson correlation
	// covariance and variances are computed up to scaling factor (n-1)
	double covXY = 0, varX = 0, varY = 0;
	for(int i = 0; i < n; i++) {
		double tmpX = x(i) - meanX, tmpY = y(i) - meanY;  // centered values
		covXY += tmpX * tmpY;
		varX += tmpX * tmpX;
		varY += tmpY * tmpY;
	}
	return covXY / (sqrt(varX) * sqrt(varY));
}


// ----------------------------------------
// Huber correlation based on winsorization
// ----------------------------------------

// C++ functions for robust correlation based on winsorization assume
// that the data are already robustly standardized

// shrink a single observation
// cm ... winsorization constant for negative values
// cp ... winsorization constant for positive values
double winsorize(const double& x, const double& cm, const double& cp) {
	if(x < cm) {
		return cm;
	} else if(x > cp) {
		return cp;
	} else return x;
}

// robust correlation based on univariate winsorization
double corHuberUni(const VectorXd& x, const VectorXd& y, const double& c) {
	const double cm = -c;		// negative winsorization constant
	const int n = x.size();
	VectorXd wx(n), wy(n);
	for(int i = 0; i < n; i++) {
		wx(i) = winsorize(x(i), cm, c);
		wy(i) = winsorize(y(i), cm, c);
	}
	// call barebones function for Pearson correlation with winsorized data
	return corPearson(wx, wy);
}

// R interface to corHuberUni()
SEXP R_corHuberUni(SEXP R_x, SEXP R_y, SEXP R_c) {
	// convert data to arma types
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	Map<VectorXd> x(Rcpp_x.begin(), Rcpp_x.size());	// reuse memory
	Map<VectorXd> y(Rcpp_y.begin(), Rcpp_y.size());	// reuse memory
	double c = as<double>(R_c);
	// call Eigen version and wrap result
	return wrap(corHuberUni(x, y, c));
}

// robust correlation based on adjusted univariate winsorization
double corHuberAdj(const VectorXd& x, const VectorXd& y, const double& c) {
	const int n = x.size();
	VectorXi sign(n);
	int n1 = 0, n2 = 0;
	for(int i = 0; i < n; i++) {
		double xy = x(i) * y(i);
		if(xy > 0) {
			sign(i) = 1;
			n1 += 1;	// count number of observations in 1. or 3. quadrant
		} else if(xy < 0) {
			sign(i) = -1;
			n2 += 1;	// count number of observations in 2. or 4. quadrant
		} else {
			sign(i) = 0;
			n1 += 1;	// borders are included
			n2 += 1;	// in both counts
		}
	}
	double c1, c2;	// constants for winsorization
	if(n1 < n2) {
		// winsorization with larger constant:  2. and 4. quadrants
		// winsorization with smaller constant: 1. and 3. quadrants
		n1 = n - n2;	// borders are included in major quadrants
		c2 = c;							// larger winsorization constant
		c1 = sqrt(n1/(double)n2) * c2;	// smaller winsorization constant
	} else {
		// winsorization with larger constant:  1. and 3. quadrants
		// winsorization with smaller constant: 2. and 4. quadrants
		// if there are the same number of observations in 1./3. and 2./4.
		// quadrants, 1./3. are assumed to be the major quadrants
		n2 = n - n1;	// borders are included in major quadrants
		c1 = c;							// larger winsorization constant
		c2 = sqrt(n2/(double)n1) * c1;	// smaller winsorization constant
	}
	// winsorize data
	const double cm = -c, cm1 = -c1, cm2 = -c2;
	VectorXd wx(n), wy(n);
	for(int i = 0; i < n; i++) {
		if(sign(i) == 1) {
			wx(i) = winsorize(x(i), cm1, c1);
			wy(i) = winsorize(y(i), cm1, c1);
		} else if(sign(i) == -1) {
			wx(i) = winsorize(x(i), cm2, c2);
			wy(i) = winsorize(y(i), cm2, c2);
		} else {
			wx(i) = winsorize(x(i), cm, c);
			wy(i) = winsorize(y(i), cm, c);
		}
	}
	// call barebones function for Pearson correlation with winsorized data
	return corPearson(wx, wy);
}

// R interface to corHuberAdj()
SEXP R_corHuberAdj(SEXP R_x, SEXP R_y, SEXP R_c) {
	// convert data to Rcpp types
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	Map<VectorXd> x(Rcpp_x.begin(), Rcpp_x.size());	// reuse memory
	Map<VectorXd> y(Rcpp_y.begin(), Rcpp_y.size());	// reuse memory
	double c = as<double>(R_c);
	// call Eigen version and wrap result
	return wrap(corHuberAdj(x, y, c));
}

// robust correlation based on bivariate winsorization
double corHuberBi(const VectorXd& x, const VectorXd& y, const double& c,
		const double& prob, const double& tol) {
	// compute the initial correlation matrix with adjusted winsorization
	double r0 = corHuberAdj(x, y, c);
	// C++ isnan() is not portable and gives error on Windows systems
	// use R macro ISNAN() instead
	if(ISNAN(r0) || (1 - abs(r0) < tol)) {
		// points almost on a line, leads to computationally singular system
		return r0;
	}
	// construct correlation matrix
	MatrixXd R0 = MatrixXd::Identity(2, 2);	// initialize identity matrix
	R0(1,0) = r0;							// set off-diagonal elements to
	R0(0,1) = r0;							// initial correlation estimate
	// compute Mahalanobis distances
	MatrixXd xy(x.size(), 2);
	xy.col(0) = x; xy.col(1) = y;	// put x and y in a matrix
	VectorXd md = (xy * R0.inverse()).cwiseProduct(xy).rowwise().sum();	// squared Mahalanobis distances
	// shrink observations with too large distances
	// compute the quantile of the chi-squared distribution with 2 degrees of
	// freedom corresponding to probability 'prob'
	// note that (1 - prob) needs to be used in the call to Rf_qchisq()
	double d = Rf_qchisq(1 - prob, 2, 0, 0);
	// shrink observations with too large distances
	for(int i = 0; i < xy.rows(); i++) {
		if(md(i) > d){
			xy.row(i) *= sqrt(d / md(i));
		}
	}
	// call barebones function for Pearson correlation with winsorized data
	return corPearson(xy.col(0), xy.col(1));
}

// R interface to corHuberBi()
SEXP R_corHuberBi(SEXP R_x, SEXP R_y, SEXP R_c, SEXP R_prob, SEXP R_tol) {
	// convert data to Rcpp types
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	Map<VectorXd> x(Rcpp_x.begin(), Rcpp_x.size());	// reuse memory
	Map<VectorXd> y(Rcpp_y.begin(), Rcpp_y.size());	// reuse memory
	double c = as<double>(R_c);
	double prob = as<double>(R_prob);
	double tol = as<double>(R_tol);
	// call Eigen version and wrap result
	return wrap(corHuberBi(x, y, c, prob, tol));
}

// robust correlation matrix based on bivariate winsorization
MatrixXd corMatHuber(const MatrixXd& x, const double& c,
		const double& prob, const double& tol) {
	const int p = x.cols();
	MatrixXd R = MatrixXd::Identity(p, p);	// initialize identity matrix
	// compute off-diagonal elements
	for(int j = 0; j < p; j++) {
		VectorXd xj = x.col(j);
		for(int i = j+1; i < p; i++) {
			R(i, j) = corHuberBi(x.col(i), xj, c, prob, tol);
			R(j, i) = R(i, j);
		}
	}
	// return robust correlation matrix
	return R;
}

// R interface to corMatHuber()
SEXP R_corMatHuber(SEXP R_x, SEXP R_c, SEXP R_prob, SEXP R_tol) {
	// convert data to Rcpp types
	NumericMatrix Rcpp_x(R_x);
	Map<MatrixXd> x(Rcpp_x.begin(), Rcpp_x.rows(), Rcpp_x.cols());	// reuse memory
	double c = as<double>(R_c);
	double prob = as<double>(R_prob);
	double tol = as<double>(R_tol);
	// call Eigen version and wrap result
	return wrap(corMatHuber(x, c, prob, tol));
}
