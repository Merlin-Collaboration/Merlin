/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef SMPTransform3D_h
#define SMPTransform3D_h 1

#include "merlin_config.h"
#include "Transform3D.h"
#include "SliceMacroParticle.h"
#include "RMap.h"

namespace SMPTracking
{

class SMPTransform3D
{
public:

	SMPTransform3D(const Transform3D& tfrm);

	/**
	 * Apply (approximate) transformation
	 */
	SliceMacroParticle& Apply(SliceMacroParticle& p) const;
	SliceMacroParticle& operator ()(SliceMacroParticle& p) const
	{
		return Apply(p);
	}

private:

	R2Map R;
	double delta_x, delta_y, theta_x, theta_y;
	bool bNoRot;
	bool nullRotation;
};

} // end namespace SMPTracking

#endif // SMPTransform3D_h
