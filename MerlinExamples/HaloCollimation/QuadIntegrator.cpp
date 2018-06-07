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

#include "TransportMatrix.h"
#include "MatrixMaps.h"
#include "PhysicalConstants.h"
#include "QuadIntegrator.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

namespace ParticleTracking
{

void QuadIntegrator::TrackStep(double ds)
{
	double dBdx = currentComponent->GetFieldStrength();
	double p0 = currentBunch->GetReferenceMomentum();
	double brho = p0 / eV / SpeedOfLight;

	RMtrx Rm(2);

	for(ParticleBunch::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
	{
		double k1 = dBdx / (brho * (1 + (*p).dp()));
		TransportMatrix::QuadrupoleR(ds, k1, Rm.R);
		Rm.Apply(*p);
	}
}

} // end namespace ParticleTracking
