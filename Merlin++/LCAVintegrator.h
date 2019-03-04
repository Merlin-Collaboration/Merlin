/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_LCAVIntegrator
#define _h_LCAVIntegrator 1

#include "TWRFStructure.h"
#include "ParticleComponentTracker.h"

namespace ParticleTracking
{

class LCAVIntegrator: public ParticleComponentTracker::Integrator<TWRFStructure>
{
protected:

	void TrackStep(double);
	void TrackEntrance();
	void TrackExit();

private:

	void ApplyEndField(double gsgn);
};

} // end namespace ParticleTracking

#endif
