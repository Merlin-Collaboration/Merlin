#include "ParticleDistributionGenerator.h"

#include "NumericalConstants.h"



PSvector NormalParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	p.x()	= RandomGauss(1,cutoffs.x());
	p.xp()	= RandomGauss(1,cutoffs.xp());
	p.y()	= RandomGauss(1,cutoffs.y());
	p.yp()	= RandomGauss(1,cutoffs.yp());
	p.dp()	= RandomGauss(1,cutoffs.dp());
	p.ct()	= RandomGauss(1,cutoffs.ct());
	return p;
}

PSvector UniformParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	p.x()	= RandomNG::uniform(-1,1);
	p.xp()	= RandomNG::uniform(-1,1);
	p.y()	= RandomNG::uniform(-1,1);
	p.yp()	= RandomNG::uniform(-1,1);
	p.dp()	= RandomNG::uniform(-1,1);
	p.ct()	= RandomNG::uniform(-1,1);
	return p;
}

PSvector RingParticleDistributionGenerator::GenerateFromDistribution() const
{
	PSvector p(0);
	double u = RandomNG::uniform(-pi,pi);
	p.x()	= cos(u);
	p.xp()	= sin(u);
	u = RandomNG::uniform(-pi,pi);
	p.y()	= cos(u);
	p.yp()	= sin(u);
	p.dp()	= RandomNG::uniform(-1,1);
	p.ct()	= RandomNG::uniform(-1,1);
	return p;
}
