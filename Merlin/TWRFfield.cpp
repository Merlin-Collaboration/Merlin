/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "TWRFfield.h"

Vector3D TWRFfield::GetBFieldAt(const Point3D& x, double t) const
{
	return Vector3D(0, 0, 0);
}

Vector3D TWRFfield::GetEFieldAt(const Point3D& x, double t) const
{
	return Vector3D(0, 0, TWRFfield::Ez(x.z, t));
}

Vector3D TWRFfield::GetForceAt(const Point3D& x, const Vector3D& v, double q, double t) const
{
	return Vector3D(0, 0, q * TWRFfield::Ez(x.z, t));
}
