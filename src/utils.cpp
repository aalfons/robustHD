/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include "utils.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// *****************
// class definitions
// *****************

class SortData {
public:
	uword index;
	double value;

	// constructors
	SortData();
	SortData(uword&, const double&);
//	// overloaded < (is less) operator for sorting and ordering
//	bool operator< (const SortData&);
};

// constructors
inline SortData::SortData() {}
inline SortData::SortData(uword& first, const double& second) {
	index = first;
	value = second;
}

//// overloaded < (is less) operator for sorting and ordering
//bool SortData::operator< (const SortData& other) {
//      return (this->value < other.value);
//}

// compare two objects with < (is less) operator for sorting and ordering
bool sortDataIsLess(const SortData& left, const SortData& right) {
	return left.value < right.value;
}


// *****************
// utility functions
// *****************

// find indices of h smallest observations
// those are not ordered, but they are the h smallest
uvec findSmallest(const vec& x, const uword& h) {
	// initialize data structure for sorting
	const uword n = x.size();
	vector<SortData> sortVector(n);
	for(uword i = 0; i < n; i++) {
		sortVector[i] = SortData(i, x(i));
	}
	// call STL's nth_element()
	nth_element(sortVector.begin(), sortVector.begin()+h, sortVector.end(),
			sortDataIsLess);
	// construct and return vector of indices
	uvec indices(h);
	for(uword i = 0; i < h; i++) {
		SortData sortElement = sortVector[i];
		indices(i) = sortElement.index;
	}
	return indices;
}

// R interface to findSmallest()
SEXP R_findSmallest(SEXP R_x, SEXP R_h) {
	NumericVector Rcpp_x(R_x);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	int h = as<int>(R_h);
	uvec indices = findSmallest(x, h);	// call native C++ function
	return wrap(indices.memptr(), indices.memptr() + indices.n_elem);
}

// find indices of h smallest observations
uvec partialOrder(const vec& x, const uword& h) {
	// initialize data structure for sorting
	const uword n = x.size();
	vector<SortData> sortVector(n);
	for(uword i = 0; i < n; i++) {
		sortVector[i] = SortData(i, x(i));
	}
	// call STL's partial_sort()
	partial_sort(sortVector.begin(), sortVector.begin()+h, sortVector.end(),
			sortDataIsLess);
	// construct and return vector of indices
	uvec indices(h);
	for(uword i = 0; i < h; i++) {
		SortData sortElement = sortVector[i];
		indices(i) = sortElement.index;
	}
	return indices;
}

// R interface to partialOrder()
SEXP R_partialOrder(SEXP R_x, SEXP R_h) {
	NumericVector Rcpp_x(R_x);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	int h = as<int>(R_h);
	uvec indices = partialOrder(x, h);	// call native C++ function
	return wrap(indices.memptr(), indices.memptr() + indices.n_elem);
}

//// R interface to partial sort of h smallest observations
//SEXP R_partialSort(SEXP R_x, SEXP R_h) {
//	NumericVector Rcpp_x(R_x);
//	NumericVector x(Rcpp_x.begin(), Rcpp_x.end());	// copy original vector
//	int h = as<int>(R_h);
//	// call STL's partial_sort()
//	partial_sort(x.begin(), x.begin()+h, x.end(), sortDataIsLess);
//	return x;
//}

// create sequence of integers (starting with 0)
uvec seqLen(const uword& n) {
	uvec sequence(n);
	for(uword i = 0; i < n; i++) {
		sequence(i) = i;
	}
	return sequence;
}

// compute sign of a numeric value
sword sign(const double& x) {
	return (x > 0) - (x < 0);
}
