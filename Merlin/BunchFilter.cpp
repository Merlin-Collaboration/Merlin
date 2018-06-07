/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "BunchFilter.h"

namespace ParticleTracking
{

ParticleBunchFilter::~ParticleBunchFilter()
{
	// Nothing to do
}

bool HorizontalHaloParticleBunchFilter::Apply(const PSvector& v) const
{
	if(v.x() > (orbit + limit) || v.x() < (orbit - limit))
	{
		if(v.x() != 0.0 && v.xp() != 0.0)
		{
			return true;
		}
	}
	return false;
}

void HorizontalHaloParticleBunchFilter::SetHorizontalLimit(double lim)
{
	cout << "Setting Horizontal limit to: " << lim << endl;
	limit = lim;
}

void HorizontalHaloParticleBunchFilter::SetHorizontalOrbit(double lim)
{
	cout << "Setting Horizontal orbit to: " << lim << endl;
	orbit = lim;
}

bool VerticalHaloParticleBunchFilter::Apply(const PSvector& v) const
{
	if(fabs(v.y()) > limit)
	{
		return true;
	}
	return false;
}

void VerticalHaloParticleBunchFilter::SetVerticalLimit(double lim)
{
	cout << "Setting Vertical limit to: " << lim << endl;
	limit = lim;
}

}
