/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_LinearAlgebra
#define _h_LinearAlgebra

#include <complex>
typedef std::complex<double> Complex;

#include "TLAS.h"

namespace TLAS
{

extern template class Vector<Complex>;
extern template class Matrix<Complex>;
extern template class SubVector<Complex>;
extern template class SubMatrix<Complex>;
extern template class ConstSubVector<Complex>;
extern template class ConstSubMatrix<Complex>;

typedef Vector<double> RealVector;
typedef Matrix<double> RealMatrix;
typedef Vector<Complex> ComplexVector;
typedef Matrix<Complex> ComplexMatrix;

/**
 * Matrix Inversion
 */
double Invert(RealMatrix& t);

/**
 * Eigensystem
 */
bool EigenSystem(RealMatrix& t, ComplexVector& eigenvalues, ComplexMatrix& eigenvectors);

/**
 * Matrix symplectification
 */
void Symplectify(RealMatrix& a);

/**
 * Eigensystem of a real symmetric matrix
 */
void EigenSystemSymmetricMatrix(RealMatrix& m, RealVector& eigenvalues);

} // end namespace TLAS;

using namespace TLAS;

#endif
