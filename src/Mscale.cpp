/*
 * Authors: Viktoria Oellerer
 *          KU Leuven
 *
 *          Andreas Alfons
 *          Erasmus Universiteit Rotterdam
 */

#include "Mscale.h"

using namespace Rcpp;
using namespace arma;


// biweight rho-function for a single value
// k ... tuning parameter
double Mrho(const double& x, const double& k) {
  if(x >= -k && x <= k) {
    return pow(x,6)/(6.0*pow(k,4)) - pow(x,4)/(2.0*pow(k,2)) + pow(x,2)/2.0;
  } else {
    return pow(k,2)/6.0;
  }
}

// biweight rho-function for a vector
// k ... tuning parameter
vec Mrho(const vec& x, const double& k) {
  const uword n = x.n_elem;
  vec out(n);
  for(uword i = 0; i < n; i++) {
    out(i) = Mrho(x(i), k);
  }
  return out;
}

// R interface to Mrho()
SEXP R_Mrho(SEXP R_x, SEXP R_k) {
  vec x = as<vec>(R_x);
  double k = as<double>(R_k);
  vec out = Mrho(x, k);
  return wrap(out.memptr(), out.memptr() + out.n_elem);
}


// biweight weight function for a single value
// k ... tuning parameter
double Mwgt(const double& x, const double& k) {
  if(x >= -k && x <= k) {
    return pow(x,4)/pow(k,4) - 2.0*pow(x,2)/pow(k,2) + 1.0;
  } else {
    return 0;
  }
}

// biweight weight function for a vector
// k ... tuning parameter
vec Mwgt(const vec& x, const double& k) {
  const uword n = x.n_elem;
  vec w(n);
  for(uword i = 0; i < n; i++) {
    w(i) = Mwgt(x(i), k);
  }
  return w;
}

// R interface to Mwgt()
SEXP R_Mwgt(SEXP R_x, SEXP R_k) {
  vec x = as<vec>(R_x);
  double k = as<double>(R_k);
  vec w = Mwgt(x, k);
  return wrap(w.memptr(), w.memptr() + w.n_elem);
}


// compute M-scale estimate
// Armadillo library is used for linear algebra
// x ....... current residuals
// k ....... tuning parameter for rho-function
// b ....... constant for consistency of M-scale
// tol ..... numerical tolerance for convergence
// nfpi .... maximum number of fixed point iterations
// warn .... logical indicating if a warning should be given if there is no
//           convergence in the maximum number of fixed point iterations
double Mscale(const vec& x, const double& k, const double& b, const uword& nfpi,
              const double& tol, const bool& warn) {
  // initializations
  uword i = 0;
  bool continueIterations = true;
  double scale = 1.4826 * median(abs(x));  // initialize scale with MAD
  // perform fixed point iterations
  while((i < nfpi) && continueIterations) {
    double previousScale = scale;
    scale *= sqrt(mean(Mrho(x/scale, k)) / b);
    double error = scale/previousScale - 1.0;
    if(error < 0) error = -error; // abs() only works for integers, not doubles
    i++;
    continueIterations = error > tol;
  }
  // give a warning  if there is no convergence in the maximum number of steps
  if(warn && (i == nfpi) && continueIterations) {
    Rf_warning("no convergence in maximum number of fixed-point iterations\n");
  }
  // return M-scale
  return scale;
}

// R interface to Mscale()
SEXP R_Mscale(SEXP R_x, SEXP R_k, SEXP R_b, SEXP R_nfpi, SEXP R_tol,
              SEXP R_warn) {
  // initializations
  vec x = as<vec>(R_x);
  double k = as<double>(R_k);
  double b = as<double>(R_b);
  double nfpi = as<double>(R_nfpi);
  double tol = as<double>(R_tol);
  bool warn = as<bool>(R_warn);
  // call native C++ function and return results
  double scale = Mscale(x, k, b, nfpi, tol, warn);
  return wrap(scale);
}
