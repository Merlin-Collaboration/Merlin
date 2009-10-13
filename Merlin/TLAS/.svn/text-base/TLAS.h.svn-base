/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:55 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_TLAS
#define _h_TLAS

// Templated functions and classes for performing linear algebra. The following
// perform linear algebra functions based on TLAS::Matrix and TLAS::Vector classes.

#include "TLAS/TMatrixLib.h"

namespace TLAS {

// Linear algebra classes.

template<class T> class LUMatrix;
template<class T> class SVDMatrix;

// Matrix inversion and determinant functions (global functions).

template<class T> Matrix<T>& InvertMatrix(Matrix<T>&);
template<class T> Matrix<T> Invert(const Matrix<T>&);
template<class T> Matrix<T> Transpose(const Matrix<T>&);
template<class T> T Det(const Matrix<T>&);

// Low level linear algebra routines (global functions).

template<class T> void ludcmp(Matrix<T>&,std::vector<int>&,T&);
template<class T, class V> V& lubksb(const Matrix<T>& a, const std::vector<int>& indx, V& b);
template<class T> void svdcmp(Matrix<T>&, Vector<T>&, Matrix<T>&);
template<class T,class V>
Vector<T>& svbksb(const Matrix<T>&, const Vector<T>&, const Matrix<T>&, const V&, Vector<T>&);

// Exceptions used by linear algebra routines.

class SingularMatrix {};
class ConvergenceFailure {};
class SingularValuesAllZero {};
class NonSquareMatrix {};


template<class T> Matrix<T> Transpose(const Matrix<T>& M)
{
    Matrix<T> Mt(M.ncols(),M.nrows());
    for(Subscript i=0; i<M.nrows(); i++)
        Mt.column(i)=M.row(i);
    return Mt;
}

}  // end namespace TLAS


#endif // _h_TLAS
