/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_RotationMatrix
#define _h_RotationMatrix 1

#include "LinearAlgebra.h"

/**
 * Constructs a 3x3 array for rotations about x
 */
inline RealMatrix DECL_XROT(double c, double s)
{
	RealMatrix m(IdentityMatrix(3));
	m(1, 1) = m(2, 2) = c;
	m(1, 2) = s;
	m(2, 1) = -s;
	return m;
}

/**
 * Constructs a 3x3 array for rotations about y
 */
inline RealMatrix DECL_YROT(double c, double s)
{
	RealMatrix m(IdentityMatrix(3));
	m(0, 0) = m(2, 2) = c;
	m(0, 2) = -s;
	m(2, 0) = s;
	return m;
}

/**
 * Constructs a 3x3 array for rotations about z
 */
inline RealMatrix DECL_ZROT(double c, double s)
{
	RealMatrix m(IdentityMatrix(3));
	m(0, 0) = m(1, 1) = c;
	m(0, 1) = s;
	m(1, 0) = -s;
	return m;
}

#endif // _h_RotationMatrix
