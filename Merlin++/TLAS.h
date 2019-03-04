/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_TLAS
#define _h_TLAS

// Templated functions and classes for performing linear algebra. The following
// perform linear algebra functions based on TLAS::Matrix and TLAS::Vector classes.

#include "TMatrixLib.h"

namespace TLAS
{

// Linear algebra classes.

template<class T> class LUMatrix;
template<class T> class SVDMatrix;

// Matrix inversion and determinant functions (global functions).

template<class T> Matrix<T>& InvertMatrix(Matrix<T>&);
template<class T> Matrix<T> Invert(const Matrix<T>&);
template<class T> Matrix<T> Transpose(const Matrix<T>&);
template<class T> T Det(const Matrix<T>&);

// Low level linear algebra routines (global functions).

template<class T> void ludcmp(Matrix<T>&, std::vector<int>&, T&);
template<class T, class V> V& lubksb(const Matrix<T>& a, const std::vector<int>& indx, V& b);
template<class T> void svdcmp(Matrix<T>&, Vector<T>&, Matrix<T>&);
template<class T, class V>
Vector<T>& svbksb(const Matrix<T>&, const Vector<T>&, const Matrix<T>&, const V&, Vector<T>&);

// Exceptions used by linear algebra routines.

class SingularMatrix
{
};
class ConvergenceFailure
{
};
class SingularValuesAllZero
{
};
class NonSquareMatrix
{
};

template<class T> Matrix<T> Transpose(const Matrix<T>& M)
{
	Matrix<T> Mt(M.ncols(), M.nrows());
	for(Subscript i = 0; i < M.nrows(); i++)
	{
		Mt.column(i) = M.row(i);
	}
	return Mt;
}

}  // end namespace TLAS

#endif // _h_TLAS
