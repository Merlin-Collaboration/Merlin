/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ParticleMapPI.h"

namespace ParticleTracking
{

void ParticleMapCI::TrackStep(double ds)
{
	assert(ds == 0);
	currentComponent->Apply(*currentBunch);
	return;
}

} // end namespace ParticleTracking
