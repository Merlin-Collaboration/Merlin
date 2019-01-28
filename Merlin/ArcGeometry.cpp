/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ArcGeometry.h"

namespace
{

inline Transform3D MakeTransform(double phi, double h)
{
	return Transform3D(Point3D((cos(phi) - 1) / h, 0, sin(phi) / h), Rotation3D::rotationY(-phi));
}

} // end anonymous namespace

Transform3D ArcGeometry::GetGeometryTransform(double s0, double s) const
{
	CheckBounds(s, s0);
	return MakeTransform((s - s0) * h, h);
}

Transform3D ArcGeometry::GetGeometryTransform(BoundaryPlane p) const
{
	double phi = p == entrance ? -GetAngle() / 2 : GetAngle() / 2;
	return MakeTransform(phi, h);
}

Transform3D ArcGeometry::GetTotalGeometryTransform() const
{
	return MakeTransform(GetAngle(), h);
}
