/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "SWRFfield.h"

Vector3D SWRFfield::GetBFieldAt(const Point3D& x, double t) const
{
	return Vector3D(0, 0, 0);
}

Vector3D SWRFfield::GetEFieldAt(const Point3D& x, double t) const
{
	return Vector3D(0, 0, SWRFfield::Ez(x.z, t));
}

Vector3D SWRFfield::GetForceAt(const Point3D& x, const Vector3D& v, double q, double t) const
{
	return Vector3D(0, 0, q * SWRFfield::Ez(x.z, t));
}
