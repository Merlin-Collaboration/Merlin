/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ElectronBunch.h"
#include "PhysicalConstants.h"

using namespace ParticleTracking;
using namespace PhysicalConstants;

double ElectronBunch::GetParticleMass() const
{
	return ElectronMass;
}

double ElectronBunch::GetParticleMassMeV() const
{
	return ElectronMassMeV;
}

double ElectronBunch::GetParticleLifetime() const
{
	return 0;
}

bool ElectronBunch::IsStable() const
{
	return true;
}
