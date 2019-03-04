/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cassert>
#include "Transformable.h"

// macro for transformations
#define _TRNSFM(func) \
	if(local_T) (*local_T) *= Transform3D::func; \
	else local_T = new Transform3D(Transform3D::func); \
	Invalidate();

Transformable::~Transformable()
{
	if(local_T)
	{
		delete local_T;
	}
}

void Transformable::Translate(double dx, double dy, double dz)
{
	_TRNSFM(translation(dx, dy, dz));
}

void Transformable::RotateX(double angle)
{
	_TRNSFM(rotationX(angle));
}

void Transformable::RotateY(double angle)
{
	_TRNSFM(rotationY(angle));
}

void Transformable::RotateZ(double angle)
{
	_TRNSFM(rotationZ(angle));
}

void Transformable::ClearTransform()
{
	if(local_T)
	{
		delete local_T;
		local_T = nullptr;
	}
	Invalidate();
}
