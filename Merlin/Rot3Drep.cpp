/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "Rot3Drep.h"
#include "IdentityRotation.h"
#include "GeneralRotation.h"
#include "AxisRotations.h"

Rot3Drep::Rot3Drep() :
	refc(0)
{
}

Rot3Drep::~Rot3Drep()
{
}

bool Rot3Drep::isIdentity() const
{
	return false;
}

bool Rot3Drep::isXrot() const
{
	return false;
}

bool Rot3Drep::isYrot() const
{
	return false;
}

bool Rot3Drep::isZrot() const
{
	return false;
}

Rot3Drep* Rot3Drep::identity()
{
	return new IdentityRotation;
}

Rot3Drep* Rot3Drep::rotationX(double angle)
{
	if(angle == 0)
	{
		return new IdentityRotation;
	}
	else
	{
		return new RotationX(angle);
	}
}

Rot3Drep* Rot3Drep::rotationY(double angle)
{
	if(angle == 0)
	{
		return new IdentityRotation;
	}
	else
	{
		return new RotationY(angle);
	}
}

Rot3Drep* Rot3Drep::rotationZ(double angle)
{
	if(angle == 0)
	{
		return new IdentityRotation;
	}
	else
	{
		return new RotationZ(angle);
	}
}
