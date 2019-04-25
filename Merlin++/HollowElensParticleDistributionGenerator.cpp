/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "HollowElensParticleDistributionGenerator.h"
#include "RandomNG.h"
#include "NumericalConstants.h"

PSvector ThickHELHaloParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi, pi);
	double randAmp = RandomNG::uniform(halo_inner, halo_outer);
	p.x()    = cos(u) * randAmp;
	p.xp()   = 0.0;
	p.y()    = sin(u) * randAmp;
	p.yp()   = 0.0;
	p.dp()   = RandomNG::uniform(-1, 1);
	p.ct()   = RandomNG::uniform(-1, 1);
	return p;
}

PSvector ThinHELHaloParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi, pi);
	p.x()    = cos(u) * halo_outer;
	p.xp()   = 0.0;
	p.y()    = sin(u) * halo_outer;
	p.yp()   = 0.0;
	p.dp()   = RandomNG::uniform(-1, 1);
	p.ct()   = RandomNG::uniform(-1, 1);
	return p;
}
