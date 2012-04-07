/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

// extend class MatrixBase with additional member functions

// append scalar to a vector
// x ... scalar to be appended
inline void append(const Scalar& x) {
	const int n = this->size();
	this->conservativeResize(n+1);
	this->operator()(n) = x;
}

// append scalar to a vector
// x ... scalar to be appended
// n ... keep n-1 observations and append scalar at position n
inline void append(const Scalar& x, const int& n) {
	this->conservativeResize(n+1);
	this->operator()(n) = x;
}


//// keep vector element
//// pos ... position of element to keep
//inline void keep(const int& pos) {
//	this->operator()(0) = this->operator()(pos);	// move element to keep
//	this->conservativeResize(1);	// resize vector
//}
//
//// keep vector elements
//// this does not work with 'Matrix<int, Dynamic, 1>', therefore it is defined
//// as template
//// pos ... positions of elements to keep (in ascending order)
//template <typename Index>
//inline void keep(const Matrix<Index, Dynamic, 1>& pos) {
//	const int k = pos.size();
//	for(int i = 0; i < k; i++) {
//		this->operator()(i) = this->operator()(pos(i));	// move elements to keep
//	}
//	this->conservativeResize(k);	// resize vector
//}
//
//// keep column of a matrix
//// pos ... positions of column to keep
//inline void keepCol(const int& pos) {
//	this->col()(0) = this->col()(pos);	// move element to keep
//	this->conservativeResize(NoChange, 1);	// resize vector
//}
//
//// keep columns of a matrix
//// this does not work with 'Matrix<int, Dynamic, 1>', therefore it is defined
//// as template
//// pos ... positions of columns to keep (in ascending order)
//template <typename Index>
//inline void keepCols(const Matrix<Index, Dynamic, 1>& pos) {
//	const int k = pos.size();
//	for(int i = 0; i < k; i++) {
//		this->col()(i) = this->col()(pos(i));	// move elements to keep
//	}
//	this->conservativeResize(NoChange, k);	// resize vector
//}

// remove vector element
// pos ... position of element to remove
inline void remove(const int& pos) {
	const int n = this->size();
	const int nKeepTail = n - (pos + 1);
	if(nKeepTail > 0) {
		// move last elements up
		this->segment(pos, nKeepTail) = this->tail(nKeepTail).eval();
	}
	this->conservativeResize(n-1);	// resize vector
}

// remove vector elements
// this does not work with 'Matrix<int, Dynamic, 1>', therefore it is defined
// as template
// pos ... positions of elements to remove (in ascending order)
template <typename Index>
inline void remove(const Matrix<Index, Dynamic, 1>& pos) {
	const int k = pos.size();
	// case of k == 0 is not considered
	if(k == 1) {
		// remove only one element
		this->remove(pos(0));
	} else {
		// remove more than one element
		const int n = this->size();
		int nFrontNew = pos(0), nFrontOld = pos(0) + 1, nKeep;
		// move intermediate blocks
		for(int i = 1; i < k; i++) {
			nKeep = pos(i) - nFrontOld;
			if(nKeep > 0) {
				this->segment(nFrontNew, nKeep) = this->segment(nFrontOld, nKeep).eval();
			}
			nFrontNew += nKeep;
			nFrontOld = pos(i) + 1;
		}
		// move last block
		nKeep = n - nFrontOld;
		if(nKeep > 0) {
			this->segment(nFrontNew, nKeep) = this->tail(nKeep).eval();
		}
		// resize vector
		this->conservativeResize(n-k);
	}
}

// remove column
// pos ... position of column to remove
inline void removeCol(const int& pos) {
//	const int n = this->rows(), p = this->cols();
	const int p = this->cols();
	const int pKeepRight  = p - (pos + 1);
	if(pKeepRight > 0) {
		// move right columns to the left
		this->middleCols(pos, pKeepRight) = this->rightCols(pKeepRight).eval();
	}
	this->conservativeResize(NoChange, p-1);	// resize matrix
}
