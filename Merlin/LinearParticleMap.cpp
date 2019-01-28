/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cassert>
#include "ParticleBunch.h"
#include "LinearParticleMap.h"

namespace ParticleTracking
{

ParticleBunch& LinearParticleMap::Apply(ParticleBunch& bunch) const
{
	R.Apply(bunch.GetParticles());
	return bunch;
}

void LinearParticleMap::Invert()
{
	// TODO:
	assert(false);
}

} //end namespace ParticleTracking
