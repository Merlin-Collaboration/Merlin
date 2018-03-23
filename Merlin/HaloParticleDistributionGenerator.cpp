#include "HaloParticleDistributionGenerator.h"
#include "RandomNG.h"
#include "NumericalConstants.h"

PSvector HorizonalHalo1ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi,pi);
	p.x()    = cos(u) * sqrt(halo_size);
	p.xp()   = sin(u) * sqrt(halo_size);
	p.y()    = 0.0;
	p.yp()   = 0.0;
	p.dp()   = RandomNG::uniform(-1,1);
	p.ct()   = RandomNG::uniform(-1,1);
	return p;
}

PSvector VerticalHalo1ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi,pi);
	p.x()    = 0.0;
	p.xp()   = 0.0;
	p.y()    = cos(u) * sqrt(halo_size);
	p.yp()   = sin(u) * sqrt(halo_size);
	p.dp()   = RandomNG::uniform(-1,1);
	p.ct()   = RandomNG::uniform(-1,1);
	return p;
}

PSvector HorizonalHalo2ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi,pi);
	p.x()    = cos(u) * sqrt(halo_size);
	p.xp()   = sin(u) * sqrt(halo_size);
	p.y()    = RandomGauss(1,cutoffs.y());
	p.yp()   = RandomGauss(1,cutoffs.yp());
	p.dp()   = RandomNG::uniform(-1,1);
	p.ct()   = RandomNG::uniform(-1,1);
	return p;
}

PSvector VerticalHalo2ParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi,pi);
	p.x()    = RandomGauss(1,cutoffs.x());
	p.xp()   = RandomGauss(1,cutoffs.xp());
	p.y()    = cos(u) * sqrt(halo_size);
	p.yp()   = sin(u) * sqrt(halo_size);
	p.dp()   = RandomNG::uniform(-1,1);
	p.ct()   = RandomNG::uniform(-1,1);
	return p;
}
