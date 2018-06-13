/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "EMField.h"

EMField::~EMField()
{
// Nothing to do
}

Vector3D EMField::GetForceAt(const Point3D& x, const Vector3D& v, double q, double t) const
{
	Vector3D B = GetBFieldAt(x, t);
	Vector3D E = GetEFieldAt(x, t);
	return q * (E + cross(v, B));
}
