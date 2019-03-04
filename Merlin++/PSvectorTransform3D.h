/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef PSvectorTransform3D_h
#define PSvectorTransform3D_h 1

#include "merlin_config.h"
#include "Transform3D.h"
#include "PSTypes.h"

/**
 *	Utility class for performing an arbitrary 3D coordinate
 *	transformation on PSvector objects. The transformation
 *	assumes that the particle(s) are fully relativistic.
 *	Since small angle approximations are assumed, large
 *	rotations about the x and y axis may lead to significant
 *	errors.
 */

class PSvectorTransform3D
{
public:

	PSvectorTransform3D(const Transform3D& tfrm);

	PSvector& Apply(PSvector& p) const;
	PSvectorArray& Apply(PSvectorArray& pv) const;
	PSvector& operator ()(PSvector& p) const;

private:

	Transform3D T;
	bool bNoRot;
};

inline PSvector& PSvectorTransform3D::operator ()(PSvector& p) const
{
	return Apply(p);
}

#endif
