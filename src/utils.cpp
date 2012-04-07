/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#include "utils.h"

using namespace Rcpp;
using namespace Eigen;
using namespace std;


// *****************
// class definitions
// *****************

class SortData {
public:
	int index;
	double value;

	// constructors
	SortData();
	SortData(int&, const double&);
//	// overloaded < (is less) operator for sorting and ordering
//	bool operator< (const SortData&);
};

// constructors
inline SortData::SortData() {}
inline SortData::SortData(int& first, const double& second) {
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
VectorXi findSmallest(const VectorXd& x, const int& h) {
	// initialize data structure for sorting
	const int n = x.size();
	vector<SortData> sortVector(n);
	for(int i = 0; i < n; i++) {
		sortVector[i] = SortData(i, x(i));
	}
	// call STL's nth_element()
	nth_element(sortVector.begin(), sortVector.begin()+h, sortVector.end(),
			sortDataIsLess);
	// construct and return vector of indices
	VectorXi indices(h);
	for(int i = 0; i < h; i++) {
		SortData sortElement = sortVector[i];
		indices(i) = sortElement.index;
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

// find indices of h smallest observations
VectorXi partialOrder(const VectorXd& x, const int& h) {
	// initialize data structure for sorting
	const int n = x.size();
	vector<SortData> sortVector(n);
	for(int i = 0; i < n; i++) {
		sortVector[i] = SortData(i, x(i));
	}
	// call STL's partial_sort()
	partial_sort(sortVector.begin(), sortVector.begin()+h, sortVector.end(),
			sortDataIsLess);
	// construct and return vector of indices
	VectorXi indices(h);
	for(int i = 0; i < h; i++) {
		SortData sortElement = sortVector[i];
		indices(i) = sortElement.index;
	}
	return indices;
}

// R interface to partialOrder()
SEXP R_partialOrder(SEXP R_x, SEXP R_h) {
	NumericVector Rcpp_x(R_x);
	Map<VectorXd> x(Rcpp_x.begin(), Rcpp_x.size());	// reuse memory
	int h = as<int>(R_h);
	VectorXi indices = partialOrder(x, h);	// call native C++ function
	return wrap(indices);
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

// order observations
VectorXi order(const VectorXd& x) {
	// initialize data structure for sorting
	const int n = x.size();
	vector<SortData> sortVector(n);
	for(int i = 0; i < n; i++) {
		sortVector[i] = SortData(i, x(i));
	}
	// call STL's sort()
	sort(sortVector.begin(), sortVector.end(), sortDataIsLess);
	// construct and return vector of indices
	VectorXi indices(n);
	for(int i = 0; i < n; i++) {
		SortData sortElement = sortVector[i];
		indices(i) = sortElement.index;
	}
	return indices;
}
