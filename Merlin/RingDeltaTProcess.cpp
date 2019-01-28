/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include <algorithm>
#include "utils.h"

#include "RingDeltaTProcess.h"
#include "SectorBend.h"

namespace ParticleTracking
{

struct ApplyDeltaT
{
private:
	double dt;

public:
	ApplyDeltaT(double _dt) :
		dt(_dt)
	{
	}

	void operator()(PSvector& v)
	{
		v.ct() += dt;
	}

};

RingDeltaTProcess::RingDeltaTProcess(int prio) :
	ParticleBunchProcess("RING DELTA T", prio)
{
}

void RingDeltaTProcess::SetCurrentComponent(AcceleratorComponent& component)
{
	active = (dynamic_cast<SectorBend*>(&component)) ? true : false;
	dL = component.GetLength();
	intS = 0;
}

void RingDeltaTProcess::DoProcess(double ds)
{
	intS += ds;
	for_each(currentBunch->begin(), currentBunch->end(), ApplyDeltaT(scale * ds));
	active = intS != dL;
}

double RingDeltaTProcess::GetMaxAllowedStepSize() const
{
	return dL - intS;
}

void RingDeltaTProcess::SetBendScale(double bendscale)
{
	scale = bendscale;
}

} // end namespace ParticleTracking
