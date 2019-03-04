/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "IdentityRotation.h"
#include "GeneralRotation.h"
#include "AxisRotations.h"

Rot3Drep* IdentityRotation::inv() const
{
	return const_cast<IdentityRotation*>(this);
}

Point3D IdentityRotation::rotate(const Point3D& p) const
{
	return p;
}

Vector3D IdentityRotation::rotate(const Vector3D& v) const
{
	return v;
}

bool IdentityRotation::isIdentity() const
{
	return true;
}

Rot3Drep* IdentityRotation::dot(const Rot3Drep& r) const
{
	return const_cast<Rot3Drep*>(&r);
}

Rot3Drep* IdentityRotation::dotBy(const RotationX& rx) const
{
	return const_cast<RotationX*>(&rx);
}

Rot3Drep* IdentityRotation::dotBy(const RotationY& ry) const
{
	return const_cast<RotationY*>(&ry);
}

Rot3Drep* IdentityRotation::dotBy(const RotationZ& rz) const
{
	return const_cast<RotationZ*>(&rz);
}

Rot3Drep* IdentityRotation::dotBy(const GeneralRotation& r) const
{
	return const_cast<GeneralRotation*>(&r);
}

RotationType IdentityRotation::type() const
{
	return ident;
}

Rot3Drep* IdentityRotation::rotXbyPI() const
{
	return const_cast<IdentityRotation*>(this);
}

Rot3Drep* IdentityRotation::rotYbyPI() const
{
	return const_cast<IdentityRotation*>(this);
}

Rot3Drep* IdentityRotation::rotZbyPI() const
{
	return const_cast<IdentityRotation*>(this);
}

RealMatrix& IdentityRotation::getMatrix(RealMatrix& m) const
{
	return m.copy(IdentityMatrix(3));
}
