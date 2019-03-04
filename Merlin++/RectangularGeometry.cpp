/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "RectangularGeometry.h"

Transform3D RectangularGeometry::GetGeometryTransform(double s0, double s) const
{
	CheckBounds(s, s0);
	return Transform3D::translation(0, 0, s - s0);
}

Transform3D RectangularGeometry::GetGeometryTransform(BoundaryPlane p) const
{
	double s = p == entrance ? -len / 2 : len / 2;
	return Transform3D::translation(0, 0, s);
}

Transform3D RectangularGeometry::GetTotalGeometryTransform() const
{
	return Transform3D::translation(0, 0, len);
}
