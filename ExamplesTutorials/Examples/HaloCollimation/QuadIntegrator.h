/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

// A quadrupole integrator which calculates the exact
// linear map for each particle energy.
//
// This local integrator is intended to override the internal
// particle integrator for multipoles, which using an second-
// order expansion for the energy.

#ifndef QuadIntegrator_h
#define QuadIntegrator_h

#include "StandardMultipoles.h"
#include "ParticleComponentTracker.h"

namespace ParticleTracking
{

class QuadIntegrator: public ParticleComponentTracker::Integrator<Quadrupole>
{
public:
	void TrackStep(double ds);
};

} // end namespace ParticleTracking

#endif
