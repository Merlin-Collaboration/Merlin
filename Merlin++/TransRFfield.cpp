/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "TransRFfield.h"

Vector3D TransverseRFfield::GetBFieldAt(const Point3D& x, double t) const
{
	return Vector3D(0, 0, 0);
}

Vector3D TransverseRFfield::GetEFieldAt(const Point3D& x, double t) const
{
	double Er = TransverseRFfield::Er(x.z, t);
	return Vector3D(Er * cos(theta), Er * sin(theta), 0);
}

Vector3D TransverseRFfield::GetForceAt(const Point3D& x, const Vector3D& v, double q, double t) const
{
	double Er = TransverseRFfield::Er(x.z, t);
	double Ex = Er * cos(theta);
	double Ey = Er * sin(theta);
	return Vector3D(q * Ex, q * Ey, 0);
}
