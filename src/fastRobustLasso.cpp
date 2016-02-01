/*
 * Authors: Andreas Alfons
 *          Erasmus Universiteit Rotterdam
 *
 *          Viktoria Oellerer
 *          KU Leuven
 */

#include "fastRobustLasso.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// **************************************************
// class definitions for sparse least trimmed squares
// **************************************************

// control parameters for sparse LTS
class SparseLTSControl {
  public:
  bool normalize;
  bool useIntercept;
  double tol;
  double eps;
  bool useGram;

  // constructor
  SparseLTSControl(const bool&, const bool&, const double&,
                   const double&, const bool&);
};

// constructor
inline SparseLTSControl::SparseLTSControl(const bool& _normalize,
                                          const bool& _useIntercept,
                                          const double& _tol,
                                          const double& _eps,
                                          const bool& _useGram) {
  normalize = _normalize;
  useIntercept = _useIntercept;
  tol = _tol;
  eps = _eps;
  useGram = _useGram;
}

// subsamples to be explored in fast-LTS algorithm
class Subsample {
  public:
  // information on the subsample
  uvec indices;
  double intercept;
  vec coefficients;
  vec residuals;
  double crit;
  bool continueSteps;

  // constructors
  Subsample();
  Subsample(const mat&, const vec&, const double&, const uvec&,
            const SparseLTSControl&);
  // compute lasso solution and residuals
  void lasso(const mat&, const vec&, const double&, const SparseLTSControl&);
  // perform C-Step
  void step(const mat&, const vec&, const double&, const SparseLTSControl&);
  // compute the residual center
  double center();
  // compute the uncorrected residual scale
  double scale(const double&);
};

// constructors
inline Subsample::Subsample() {
  crit = R_PosInf;
  continueSteps = true;
}
inline Subsample::Subsample(const mat& x, const vec& y,
                            const double& lambda,
                            const uvec& initial,
                            const SparseLTSControl& control) {
  const uword h = initial.n_elem;
  indices = uvec(h);
  for(uword i = 0; i < h; i++) {
    indices(i) = initial(i);	// data are copied
  }
  lasso(x, y, lambda, control); // compute lasso fit on initial subsample
  continueSteps = true;
}

// compute lasso solution, residuals and value of objective function
void Subsample::lasso(const mat& x, const vec& y, const double& lambda,
                      const SparseLTSControl& control) {
  // call standalone function
  fastLasso(x, y, lambda, true, indices, control.normalize,
            control.useIntercept, control.eps, control.useGram,
            true, intercept, coefficients, residuals, crit);
}

// perform C-Step
void Subsample::step(const mat& x, const vec& y, const double& lambda,
                     const SparseLTSControl& control) {
  // update subset
  const uword h = indices.n_elem;
  indices = findSmallest(abs(residuals), h);
  // compute lasso solution for new subset
  // (also updated residuals and value of objective function)
  double previousCrit = crit;
  lasso(x, y, lambda, control);
  // check for convergence
  continueSteps = ((previousCrit - crit) > control.tol);
}

// compute the residual center
double Subsample::center() {
  const uword h = indices.n_elem;
  double mean = 0;
  for(uword i = 0; i < h; i++) {
    mean += residuals(indices(i));
  }
  mean /= h;
  return mean;
}

// compute the uncorrected residual scale on the h smallest observations
// (center estimate is passed as parameter)
double Subsample::scale(const double& center) {
  // initialize STL vector for sorting
  const int n = residuals.n_elem, h = indices.n_elem;
  vector<double> squares(n);
  for(int i = 0; i < n; i++) {
    squares[i] = pow(residuals(i)-center, 2); // squared centered residuals
  }
  // call STL's nth_element()
  nth_element(squares.begin(), squares.begin()+h, squares.end());
  // compute scale estimate of h smallest residuals
  double sumOfSquares = 0;
  for(int i = 0; i < h; i++) {
    sumOfSquares += squares[i];
  }
  return sqrt(sumOfSquares / (double)h);
}

// compare two objects with < (is less) operator for sorting and ordering
bool subsampleIsLess(const Subsample& left, const Subsample& right) {
  return left.crit < right.crit;
}

// compare two objects with == (is equal) operator
bool subsampleIsEqual(const Subsample& left, const Subsample& right) {
  bool equal = (left.crit == right.crit);
  if(equal) {
    // values of the objective function are equal
    // check if indices are equal too
    uvec leftIndices = sort(left.indices);
    uvec rightIndices = sort(right.indices);
    // loop over sorted indices to see if they are all equal
    uword i = 0, h = leftIndices.n_elem;
    while(equal && (i < h)) {
      equal = (leftIndices(i) == rightIndices(i));
      i++;
    }
  }
  return equal;
}

// keep best subsamples
void keepBest(vector<Subsample>& subsamples, int& nkeep) {
  // sort subsamples
  sort(subsamples.begin(), subsamples.end(), subsampleIsLess);
  // loop over subsamples until nkeep unique subsamples are found
  // k counts the number of unique subsamples found
  int k = 1, n = subsamples.size();
  while((k < nkeep) && (k < n)) {
    if(subsampleIsEqual(subsamples[k-1], subsamples[k])) {
      subsamples.erase(subsamples.begin()+k);
      n--;
    } else {
      k++;
    }
  }
  // adjust nkeep if there are not enough unique subsamples
  if(k < nkeep) {
    nkeep = k;
  }
  // resize vector to keep the best subsamples
  subsamples.resize(nkeep);
}


// ****************************
// sparse least trimmed squares
// ****************************

// R interface for raw lasso fit on full subset
// (used if R function sparseLTS() is called without trimming)
SEXP R_rawLasso(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_normalize,
                SEXP R_intercept, SEXP R_eps, SEXP R_useGram) {
  // data initializations
  NumericMatrix Rcpp_x(R_x);        	// predictor matrix
  const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
  mat x(Rcpp_x.begin(), n, p, false);	// reuse memory
  NumericVector Rcpp_y(R_y);			    // response
  vec y(Rcpp_y.begin(), n, false);	  // reuse memory
  double lambda = as<double>(R_lambda);
  bool normalize = as<bool>(R_normalize);
  bool useIntercept = as<bool>(R_intercept);
  double tol = R_NaReal; // not used, but needs to be defined for control object
  double eps = as<double>(R_eps);
  bool useGram = as<bool>(R_useGram);
  SparseLTSControl control(normalize, useIntercept, tol, eps, useGram);
  // indices of full data set
  uvec indices(n);
  for(int i = 0; i < n; i++) indices(i) = i;
  // initialize object for subset, which calls native C++ function to fit lasso
  Subsample subsample(x, y, lambda, indices, control);
  double center = subsample.center();      // residual center
  double scale = subsample.scale(center);  // uncorrected residual scale
  // return results as list
  vec coefficients = subsample.coefficients;
  if(useIntercept) {
    // prepend intercept
    coefficients.insert_rows(0, 1, false);
    coefficients(0) = subsample.intercept;
  }
  return List::create(
    Named("best") = subsample.indices + 1,
    Named("coefficients") = coefficients,
    Named("residuals") = subsample.residuals,
    Named("objective") = subsample.crit,
    Named("center") = center,
    Named("scale") = scale
    );
}

// sparse least trimmed squares
// Armadillo library is used for linear algebra
// x ......... predictor matrix
// y ......... response
// lambda .... penalty parameter
// initial ... matrix of initial subsets
// nstep ..... number of initial C-steps
// nkeep ..... number of subsets to keep after initial C-steps
// control ... control parameters for lasso fits on subsamples
// ncores .... number of processor cores for parallel computing
Subsample fastSparseLTS(const mat& x, const vec& y, const double& lambda,
                        const umat& initial, const int& nstep, int& nkeep,
                        const SparseLTSControl& control, int& ncores) {

  // initializations
  const int nsamp = initial.n_cols;

  // block for parallel computing
  vector<Subsample> resamples(nsamp);
  #pragma omp parallel num_threads(ncores)
  {
    // define STL vector of initial subsamples, compute lasso fits and
    // perform initial C-steps
    #pragma omp for schedule(dynamic)
    for(int k = 0; k < nsamp; k++) {
      // compute lasso fit on initial subsample
      Subsample resampleK(x, y, lambda, initial.unsafe_col(k), control);
      // perform initial C-steps
      int i = 0;
      while(resampleK.continueSteps && (i < nstep)) {
        resampleK.step(x, y, lambda, control);
        i++;
      }
      resamples[k] = resampleK;
    }

    // keep best subsamples (let only one process do the work)
    #pragma omp single
    if(nkeep < nsamp) {
      keepBest(resamples, nkeep);
    }

    // perform additional C-steps on best subsamples until convergence
    #pragma omp for schedule(dynamic)
    for(int k = 0; k < nkeep; k++) {
      Subsample resampleK = resamples[k];
      while(resampleK.continueSteps) {
        resampleK.step(x, y, lambda, control);
      }
      resamples[k] = resampleK;
    }
  }

  // find optimal subsample
  int which = 0;              // initial optimal subsample
  double minCrit = R_PosInf;  // initial minimum of objective function
  for(int k = 0; k < nkeep; k++) {
    Subsample resampleK = resamples[k];
    if(resampleK.crit < minCrit) {
      which = k;                // update optimal subsample
      minCrit = resampleK.crit; // update minimum of objective function
    }
  }

  // keep best resample
  Subsample best = resamples[which];
  return best;
}

// R interface to fastSparseLTS()
// initial subsamples are constructed in R and passed down to C++
SEXP R_fastSparseLTS(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_initial,
                     SEXP R_normalize, SEXP R_intercept, SEXP R_ncstep,
                     SEXP R_nkeep, SEXP R_tol, SEXP R_eps, SEXP R_useGram,
                     SEXP R_ncores) {

  // data initializations
  NumericMatrix Rcpp_x(R_x);					// predictor matrix
  const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
  mat x(Rcpp_x.begin(), n, p, false);	// reuse memory
  NumericVector Rcpp_y(R_y);			    // response
  vec y(Rcpp_y.begin(), n, false);	  // reuse memory
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
  SparseLTSControl control(normalize, useIntercept, tol, eps, useGram);
  int ncores = as<int>(R_ncores);
  // call native C++ function
  Subsample best = fastSparseLTS(x, y, lambda, initial, ncstep, nkeep,
                                 control, ncores);
  double center = best.center();      // residual center
  double scale = best.scale(center);  // uncorrected residual scale

  // prepend intercept if necessary
  vec coefficients = best.coefficients;
  if(useIntercept) {
    coefficients.insert_rows(0, 1, false);
    coefficients(0) = best.intercept;
  }

  // return results as list
  return List::create(
    Named("best") = best.indices + 1,
    Named("coefficients") = coefficients,
    Named("residuals") = best.residuals,
    Named("objective") = best.crit,
    Named("center") = center,
    Named("scale") = scale
    );
}


// ****************************************
// class definitions for sparse S-estimator
// ****************************************

// control parameters for sparse S
class SparseSControl {
  public:
  // control parameters for M-estimator of scale
  double k;
  double b;
  uword nfpi;
  uword nfpiMax;
  // control parameters for weighted lasso
  bool findLambda;
  double tol;
  double eps;
  bool useGram;

  // constructors
  SparseSControl(const double&, const double&, const uword&, const uword&,
                 const bool&, const double&, const double&, const bool&);
};

// constructor
inline SparseSControl::SparseSControl(const double& _k,
                                      const double& _b,
                                      const uword& _nfpi,
                                      const uword& _nfpiMax,
                                      const bool& _findLambda,
                                      const double& _tol,
                                      const double& _eps,
                                      const bool& _useGram) {
  k = _k;
  b = _b;
  nfpi = _nfpi;
  nfpiMax = _nfpiMax;
  findLambda = _findLambda;
  tol = _tol;
  eps = _eps;
  useGram = _useGram;
}

// biweight weight function for sparse S-estimator
// k ... tuning parameter
vec Swgt(const vec& x, const double& k) {
  vec w0 = Mwgt(x, k);              // weights of unpenalized S-estimator
  return w0 / mean(w0 % pow(x, 2)); // weights of sparse S-estimator
}

// weighted samples to be explored in fast-S algorithm
class WeightedSample {
  public:
  // information on the weighted sample
  double lambda;
  vec weights;
  double intercept;
  vec coefficients;
  vec residuals;
  double scale;
  double crit;
  bool continueSteps;

  // constructors
  WeightedSample();
  WeightedSample(const mat&, const vec&, const double&, const uvec&,
                 const SparseSControl&);
  // compute lasso solution and residuals
  void lasso(const mat&, const vec&, const SparseSControl&);
  // perform I-Step
  void step1(const mat&, const vec&, const SparseSControl&);
  void step2(const mat&, const vec&, const SparseSControl&);
};

// constructors
inline WeightedSample::WeightedSample() {
  // if the optimal lambda should be found in the weighted lasso fits, it is
  // important to initialize it with infinity.  this way, when computing the
  // initial lasso with three observations, the most important variable always
  // enters the fit.
  lambda = R_PosInf;
  crit = R_PosInf;
  continueSteps = true;
}
inline WeightedSample::WeightedSample(const mat& x, const vec& y,
                                      const double& fixedLambda,
                                      const uvec& initial,
                                      const SparseSControl& control) {
  if(!control.findLambda) lambda = fixedLambda;
  // compute lasso fit on initial subsample with three observations
  fastLasso(x, y, lambda, true, initial, false, true, control.eps,
            control.useGram, false, intercept, coefficients, residuals, crit);
  // compute initial scale and weights
  scale = Mscale(residuals, control.k, control.b, control.nfpi,
                 control.tol, false);
  weights = Swgt(residuals / scale, control.k);
  // compute weighted lasso
  lasso(x, y, control);
  // compute residual M-scale with only a few fixed point iterations
  scale = Mscale(residuals, control.k, control.b, control.nfpi,
                 control.tol, false);
  // compute value of objective function
  crit = pow(scale, 2) + lambda * norm(coefficients, 1);
  continueSteps = true;
}

// compute lasso solution and residuals
void WeightedSample::lasso(const mat& x, const vec& y,
                           const SparseSControl& control) {
  // call standalone function
  fastLasso(x, y, control.findLambda, lambda, true, weights, control.tol,
            control.eps, control.useGram, intercept, coefficients, residuals);
}

// perform I-Step
// limited number of fixed-point iterations in first phase
void WeightedSample::step1(const mat& x, const vec& y,
                          const SparseSControl& control) {
  // update weights
  weights = Swgt(residuals / scale, control.k);
  // compute lasso solution for new weights (also updated residuals)
  lasso(x, y, control);
  // update scale
  scale = Mscale(residuals, control.k, control.b, control.nfpi,
                 control.tol, false);
  // compute value of objective function
  double previousCrit = crit;
  crit = pow(scale, 2) + lambda * norm(coefficients, 1);
  // check for convergence
  continueSteps = ((previousCrit - crit) > control.tol);
}
// fully iterated M-scales in second phase
void WeightedSample::step2(const mat& x, const vec& y,
                           const SparseSControl& control) {
  // update weights
  weights = Swgt(residuals / scale, control.k);
  // compute lasso solution for new weights (also updated residuals)
  lasso(x, y, control);
  // update scale
  scale = Mscale(residuals, control.k, control.b, control.nfpiMax,
                 control.tol, false);
  // compute value of objective function
  double previousCrit = crit;
  crit = pow(scale, 2) + lambda * norm(coefficients, 1);
  // check for convergence
  continueSteps = ((previousCrit - crit) > control.tol);
}

// compare two objects with < (is less) operator for sorting and ordering
bool weightedSampleIsLess(const WeightedSample& left,
                          const WeightedSample& right) {
  return left.crit < right.crit;
}

// keep best subsamples
void keepBest(vector<WeightedSample>& weightedSamples, const int& nkeep) {
  // find the best weighted samples and put them first in the vector
  nth_element(weightedSamples.begin(), weightedSamples.begin()+nkeep,
              weightedSamples.end(), weightedSampleIsLess);
  // resize vector to keep the best subsamples
  weightedSamples.resize(nkeep);
}


// ******************
// sparse S-estimator
// ******************

// TODO: properly share code between sparse LTS and sparse S
// This should be done with a shared workhorse function for the fast-LTS and
// fast-S algorithm.  However, this is complicated because M-scales are not
// fully iterated in the first phase of the fast-S algorithm.  Hence the first
// phase should always be completed for all resamples (even if there is
// convergence before) such that the M-scale is fully iterated.  In addition,
// the control object for sparse S needs to be updated in the last iteration of
// the first phase such that the M-scales are fully iterated from then on.

// sparse S-estimator
// Armadillo library is used for linear algebra
// x ............. predictor matrix
// y ............. response
// fixedLambda ... penalty parameter
// initial ....... matrix of initial subsets
// nstep ......... number of initial I-steps
// nkeep ......... number of subsets to keep after initial I-steps
// control ....... control parameters for weighted lasso and M-scale
// ncores ........ number of processor cores for parallel computing
WeightedSample fastSparseS(const mat& x, const vec& y,
                           const double& fixedLambda, const umat& initial,
                           const int& nstep, const int& nkeep,
                           SparseSControl& control, int& ncores) {

  // initializations
  const int nsamp = initial.n_cols;

  // block for parallel computing
  vector<WeightedSample> resamples(nsamp);
  #pragma omp parallel num_threads(ncores)
  {
    // define STL vector of initial resamples, compute lasso fits and
    // perform initial I-steps
    #pragma omp for schedule(dynamic)
    for(int k = 0; k < nsamp; k++) {
      // compute lasso fit on initial resample
      WeightedSample resampleK(x, y, fixedLambda, initial.unsafe_col(k),
                               control);
      // perform initial I-steps with only a few iterations for M-scale
      // last initial I-step requires fully iterated M-scale
      int i = 0, nloop = nstep-1;
      while(i < nloop) {
        resampleK.step1(x, y, control);
        i++;
      }
      // last initial I-step with fully iterated M-scale
      resampleK.step2(x, y, control);
      // save updated resample in vector
      resamples[k] = resampleK;
    }

    // keep best resamples (let only one process do the work)
    #pragma omp single
    if(nkeep < nsamp) {
      keepBest(resamples, nkeep);
    }

    // perform additional I-steps on best subsets until convergence
    #pragma omp for schedule(dynamic)
    for(int k = 0; k < nkeep; k++) {
      WeightedSample resampleK = resamples[k];
      while(resampleK.continueSteps) {
        resampleK.step2(x, y, control);
      }
      resamples[k] = resampleK;
    }
  }

  // find optimal resample
  int which = 0;              // initial optimal resample
  double minCrit = R_PosInf;  // initial minimum of objective function
  for(int k = 0; k < nkeep; k++) {
    WeightedSample resampleK = resamples[k];
    if(resampleK.crit < minCrit) {
      which = k;                // update optimal resample
      minCrit = resampleK.crit; // update minimum of objective function
    }
  }

  // keep best resample
  WeightedSample best = resamples[which];
  return best;
}

// R interface to fastSparseS()
// initial subsets are constructed in R and passed down to C++
SEXP R_fastSparseS(SEXP R_x, SEXP R_y, SEXP R_findLambda, SEXP R_fixedLambda,
                   SEXP R_initial, SEXP R_nistep, SEXP R_nkeep, SEXP R_k,
                   SEXP R_b, SEXP R_nfpi, SEXP R_nfpiMax, SEXP R_tol,
                   SEXP R_eps, SEXP R_useGram, SEXP R_ncores) {

  // data initializations
  NumericMatrix Rcpp_x(R_x);  				// predictor matrix
  const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
  mat x(Rcpp_x.begin(), n, p, false);	// reuse memory
  NumericVector Rcpp_y(R_y);			    // response
  vec y(Rcpp_y.begin(), n, false);	  // reuse memory
  bool findLambda = as<bool>(R_findLambda);
  double fixedLambda = as<double>(R_fixedLambda);
  IntegerMatrix Rcpp_initial(R_initial);	// matrix of initial subsets
  const int h = Rcpp_initial.nrow(), nsamp = Rcpp_initial.ncol();
  umat initial(h, nsamp);
  for(int j = 0; j < nsamp; j++) {
    for(int i = 0; i < h; i++) {
      // can't use the same memory-saving conversion for integer matrices
      initial(i,j) = Rcpp_initial(i,j) - 1;
    }
  }
  int nistep = as<int>(R_nistep);
  int nkeep = as<int>(R_nkeep);
  double k = as<double>(R_k);
  double b = as<double>(R_b);
  uword nfpi = as<uword>(R_nfpi);
  uword nfpiMax = as<uword>(R_nfpiMax);
  double tol = as<double>(R_tol);
  double eps = as<double>(R_eps);
  bool useGram = as<bool>(R_useGram);
  SparseSControl control(k, b, nfpi, nfpiMax, findLambda, tol, eps, useGram);
  int ncores = as<int>(R_ncores);
  // call native C++ function
  WeightedSample best = fastSparseS(x, y, fixedLambda, initial, nistep, nkeep,
                                    control, ncores);

  // prepend intercept
  vec coefficients = best.coefficients;
  coefficients.insert_rows(0, 1, false);
  coefficients(0) = best.intercept;

  // return results as list
  return List::create(
    Named("lambda") = best.lambda,
    Named("weights") = best.weights,
    Named("coefficients") = coefficients,
    Named("residuals") = best.residuals,
    Named("scale") = best.scale,
    Named("objective") = best.crit
    );
}


// *******************************************
// sparse MM-estimator
// *******************************************
// computation of sparse MM-estimate
// Armadillo library is used for linear algebra
// x .............. predictor matrix
// y .............. response
// k .............. tuning parameter
// rmax............ maximum number of iteration in weighted lasso estimation
// control......... object including parameters
// best............ object including initial fit and is in the end used to return solution

void fastSparseMM(const mat& x, const vec& y, const double& fixedLambda,
                  const double fixedLambdaS, const double& k, const int& rMax,
                  const umat& initial, const int& nstep, const int& nkeep,
                  SparseSControl& control, int& ncores,
                  // solutions of sparse S and sparse MM are returned through
                  // the following parameters
                  WeightedSample& best, WeightedSample& bestS) {

  // compute initial S-estimator
  bestS = fastSparseS(x, y, fixedLambdaS, initial, nstep, nkeep,
                      control, ncores);

  // initialize solution of MM-estimator
  best = bestS;
  best.scale = mean(Mrho(best.residuals/bestS.scale, k)) * pow(bestS.scale, 2);
  // another criterion is used for sparse MM than for sparse S
  if(control.findLambda) {
    // cannot use the lambda found in S-estimator since penalty parameter plays
    // a different role
    // TODO: call R function lambdaToLambdaS()
    best.lambda = R_PosInf;
    best.crit = R_PosInf;
  } else {
    best.lambda = fixedLambda;
    best.crit = best.scale + best.lambda * norm(best.coefficients, 1);
  }
  best.continueSteps = true;

  // reweighting steps (iteratively reweighted lasso)
  int r = 0;
  double previousCrit;
  while((r < rMax) && best.continueSteps) {
    previousCrit = best.crit;
    // update weights
    best.weights = Mwgt(best.residuals / bestS.scale, k);
    best.lambda *= 2.0; // first order condition requires 2*lambda
    fastLasso(x, y, control.findLambda, best.lambda, true, best.weights,
              control.tol, control.eps, control.useGram, best.intercept,
              best.coefficients, best.residuals);
    best.lambda /= 2.0; // transform lambda back to actual value
    // compute value of objective function
    best.scale = mean(Mrho(best.residuals/bestS.scale, k)) * pow(bestS.scale,2);
    best.crit = best.scale + best.lambda * norm(best.coefficients, 1);
    // check for convergence
    best.continueSteps = (abs(previousCrit - best.crit) > control.tol);
    r++;
  }

  // warning if there is no convergence in maximum number of iteration steps
  if((r == rMax) && best.continueSteps) {
    Rf_warning("no convergence in maximum number of reweighting steps\n");
  }
}

// R interface to fastSparseMM()
// initial subsets are constructed in R and passed down to C++
SEXP R_fastSparseMM(SEXP R_x, SEXP R_y, SEXP R_findLambda, SEXP R_fixedLambda,
                    SEXP R_fixedLambdaS, SEXP R_k, SEXP R_rMax, SEXP R_initial,
                    SEXP R_nistep, SEXP R_nkeep, SEXP R_kS, SEXP R_b,
                    SEXP R_nfpi, SEXP R_nfpiMax, SEXP R_tol, SEXP R_eps,
                    SEXP R_useGram, SEXP R_ncores) {

  // data initializations
  NumericMatrix Rcpp_x(R_x);    			// predictor matrix
  const int n = Rcpp_x.nrow(), p = Rcpp_x.ncol();
  mat x(Rcpp_x.begin(), n, p, false);	// reuse memory
  NumericVector Rcpp_y(R_y);			    // response
  vec y(Rcpp_y.begin(), n, false);	  // reuse memory
  bool findLambda = as<bool>(R_findLambda);
  double fixedLambda = as<double>(R_fixedLambda);
  double fixedLambdaS = as<double>(R_fixedLambdaS);
  double k = as<double>(R_k);
  int rMax = as<int>(R_rMax);
  IntegerMatrix Rcpp_initial(R_initial);	// matrix of initial subsets
  const int h = Rcpp_initial.nrow(), nsamp = Rcpp_initial.ncol();
  umat initial(h, nsamp);
  for(int j = 0; j < nsamp; j++) {
    for(int i = 0; i < h; i++) {
      // can't use the same memory-saving conversion for integer matrices
      initial(i,j) = Rcpp_initial(i,j) - 1;
    }
  }
  int nistep = as<int>(R_nistep);
  int nkeep = as<int>(R_nkeep);
  double kS = as<double>(R_kS);
  double b = as<double>(R_b);
  uword nfpi = as<uword>(R_nfpi);
  uword nfpiMax = as<uword>(R_nfpiMax);
  double tol = as<double>(R_tol);
  double eps = as<double>(R_eps);
  bool useGram = as<bool>(R_useGram);
  SparseSControl control(kS, b, nfpi, nfpiMax, findLambda, tol, eps, useGram);
  int ncores = as<int>(R_ncores);
  // call native C++ function
  WeightedSample best, bestS;
  fastSparseMM(x, y, fixedLambda, fixedLambdaS, k, rMax, initial, nistep,
               nkeep, control, ncores, best, bestS);

  // prepend intercept
  // sparse MM-estimator
  vec coefficients = best.coefficients;
  coefficients.insert_rows(0, 1, false);
  coefficients(0) = best.intercept;
  // sparse S-estimator
  vec coefficientsS = bestS.coefficients;
  coefficientsS.insert_rows(0, 1, false);
  coefficientsS(0) = bestS.intercept;

  // return results as list
  return List::create(
    Named("lambda") = best.lambda,
    Named("weights") = best.weights,
    Named("coefficients") = coefficients,
    Named("residuals") = best.residuals,
    Named("scale") = best.scale,
    Named("objective") = best.crit,
    Named("lambdaS") = bestS.lambda,
    Named("weightsS") = bestS.weights,
    Named("coefficientsS") = coefficientsS,
    Named("residualsS") = bestS.residuals,
    Named("scaleS") = bestS.scale,
    Named("objectiveS") = bestS.crit
    );
}
