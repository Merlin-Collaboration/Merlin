/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "HaloParticleDistributionGenerator.h"
#include "RandomNG.h"
#include "NumericalConstants.h"

PSvector HorizonalHalo1ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi, pi);
	p.x()    = cos(u) * halo_size;
	p.xp()   = sin(u) * halo_size;
	p.y()    = 0.0;
	p.yp()   = 0.0;
	p.dp()   = RandomNG::uniform(-1, 1);
	p.ct()   = RandomNG::uniform(-1, 1);
	return p;
}

PSvector VerticalHalo1ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi, pi);
	p.x()    = 0.0;
	p.xp()   = 0.0;
	p.y()    = cos(u) * halo_size;
	p.yp()   = sin(u) * halo_size;
	p.dp()   = RandomNG::uniform(-1, 1);
	p.ct()   = RandomNG::uniform(-1, 1);
	return p;
}

PSvector HorizonalHalo2ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi, pi);
	p.x()    = cos(u) * halo_size;
	p.xp()   = sin(u) * halo_size;
	p.y()    = RandomGauss(1, cutoffs.y());
	p.yp()   = RandomGauss(1, cutoffs.yp());
	p.dp()   = RandomNG::uniform(-1, 1);
	p.ct()   = RandomNG::uniform(-1, 1);
	return p;
}

PSvector VerticalHalo2ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi, pi);
	p.x()    = RandomGauss(1, cutoffs.x());
	p.xp()   = RandomGauss(1, cutoffs.xp());
	p.y()    = cos(u) * halo_size;
	p.yp()   = sin(u) * halo_size;
	p.dp()   = RandomNG::uniform(-1, 1);
	p.ct()   = RandomNG::uniform(-1, 1);
	return p;
}
