/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include "MuonBunch.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "MerlinException.h"
#include "Aperture.h"

using namespace std;
using namespace ParticleTracking;
using namespace PhysicalConstants;
using namespace PhysicalUnits;

double MuonBunch::GetParticleMass() const
{
	return MuonMass;
}

double MuonBunch::GetParticleMassMeV() const
{
	return MuonMassMeV;
}

double MuonBunch::GetParticleLifetime() const
{
	return MuonLifetime;
}

bool MuonBunch::IsStable() const
{
	return false;
}
