/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "TMatrixLib.h"

#include <complex>
typedef std::complex<double> Complex;

namespace TLAS
{
template class Vector<double>;
template class Matrix<double>;
template class SubVector<double>;
template class SubMatrix<double>;
template class ConstSubVector<double>;
template class ConstSubMatrix<double>;
}
