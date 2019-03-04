/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ConstantStrayFieldProcess.h"
#include "AcceleratorComponent.h"
#include "PhysicalConstants.h"
#include "PhysicalUnits.h"
#include "utils.h"
#include <algorithm>

using namespace PhysicalConstants;
using namespace PhysicalUnits;

namespace ParticleTracking
{

// construction
ConstantStrayFieldProcess::ConstantStrayFieldProcess(double mxstp, double b) :
	ParticleBunchProcess("Constant Stray Field Proc", 99), By(b), maxStep(mxstp), kick_ds(0), s_int(0), cL(0)
{
}

void ConstantStrayFieldProcess::SetCurrentComponent(AcceleratorComponent& component)
{
	// We override this function so that we know when we start to track a new component.
	// We are not interested in the type of component, as we apply the kick irrespective of
	// the type.

	// First call parent method
	ParticleBunchProcess::SetCurrentComponent(component);

	cL = component.GetLength();
	kick_ds = _MIN(maxStep, cL);
	s_int = 0;
	active =  true; // this process is always active.
}

void ConstantStrayFieldProcess::DoProcess(double ds)
{
	s_int += ds;
	cL -= ds;
	if(fequal(s_int, kick_ds))
	{
		ApplyKick();
		// calculate next kick_ds
		kick_ds = _MIN(maxStep, cL);
		s_int = 0;
	}
}

double ConstantStrayFieldProcess::GetMaxAllowedStepSize() const
{
	return kick_ds - s_int;
}

void ConstantStrayFieldProcess::ApplyKick()
{
	double brho = (currentBunch->GetReferenceMomentum()) / eV / SpeedOfLight;
	double dxp  = -(currentBunch->GetChargeSign()) * By * kick_ds / brho;
	for(ParticleBunch::iterator p = currentBunch->begin(); p != currentBunch->end(); p++)
	{
		p->xp() += dxp;
	}
}

} // end namespace ParticleTracking
