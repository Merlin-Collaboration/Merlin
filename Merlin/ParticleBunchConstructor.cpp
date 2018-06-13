/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <memory>

#include "RandomNG.h"
#include "ParticleBunchConstructor.h"
#include "NormalTransform.h"
#include "NumericalConstants.h"

namespace ParticleTracking
{

inline double RandomGauss(double variance, double cutoff)
{
	return cutoff == 0 ? RandomNG::normal(0, variance) :  RandomNG::normal(0, variance, cutoff);
}

ParticleBunchConstructor::ParticleBunchConstructor(const BeamData& beam, size_t npart, DistributionType dist) :
	np(npart), dtype(dist), cutoffs(0), beamdat(beam), itsFilter(nullptr), M(NormalTransform(beam)), force_c(false)
{
}

ParticleBunchConstructor::~ParticleBunchConstructor()
{
	if(itsFilter)
	{
		delete itsFilter;
	}
}

void ParticleBunchConstructor::SetBunchData(const BeamData& beam)
{
	beamdat = beam;
	M.R = NormalTransform(beam);
}

void ParticleBunchConstructor::SetNumParticles(size_t npart)
{
	assert(npart > 0);
	np = npart;
}

void ParticleBunchConstructor::SetDistributionCutoff(double cut)
{
	cutoffs = PSvector(fabs(cut));
}

void ParticleBunchConstructor::SetDistributionCutoff(const PSvector& cut)
{
	cutoffs = cut;
}

void ParticleBunchConstructor::ConstructBunchDistribution(int bunchIndex) const
{
	// The first particle is *always* the centroid particle
	PSvector p;
	p.x() = beamdat.x0;
	p.xp() = beamdat.xp0;
	p.y() = beamdat.y0;
	p.yp() = beamdat.yp0;
	p.dp() = 0;
	p.ct() = beamdat.ct0;
	p.type() = -1.0;
	p.location() = -1.0;
	p.id() = 0;
	p.sd() = 0.0;
	pbunch.push_back(p);

	size_t i = 1;
	while(i < np)
	{
		p = GenerateFromDistribution();

		// apply emittance
		p.x() *= sqrt(beamdat.emit_x);
		p.xp() *= sqrt(beamdat.emit_x);
		p.y() *= sqrt(beamdat.emit_y);
		p.yp() *= sqrt(beamdat.emit_y);
		p.dp() *= sqrt(beamdat.sig_dp);
		p.ct() *= sqrt(beamdat.sig_z);

		// Apply Courant-Snyder
		M.Apply(p);

		p += pbunch.front(); // add centroid

		p.type() = -1.0;
		p.location() = -1.0;
		p.id() = i;
		p.sd() = 0.0;

		if(itsFilter == nullptr || itsFilter->Apply(p))
		{
			pbunch.push_back(p);
			i++;
		}
	}

	if(force_c)
	{
		DoForceCentroid();
	}
}

PSvector ParticleBunchConstructor::GenerateFromDistribution() const
{
	PSvector p;
	double u;
	switch(dtype)
	{
	case normalDistribution:
		p.x()   = RandomGauss(1, cutoffs.x());
		p.xp()  = RandomGauss(1, cutoffs.xp());
		p.y()   = RandomGauss(1, cutoffs.y());
		p.yp()  = RandomGauss(1, cutoffs.yp());
		p.dp()  = RandomGauss(1, cutoffs.dp());
		p.ct()  = RandomGauss(1, cutoffs.ct());
		break;
	case flatDistribution:
		p.x()   = RandomNG::uniform(-1, 1);
		p.xp()  = RandomNG::uniform(-1, 1);
		p.y()   = RandomNG::uniform(-1, 1);
		p.yp()  = RandomNG::uniform(-1, 1);
		p.dp()  = RandomNG::uniform(-1, 1);
		p.ct()  = RandomNG::uniform(-1, 1);
		break;
	case skewHaloDistribution:
	case ringDistribution:
		u = RandomNG::uniform(-pi, pi);
		p.x()   = cos(u);
		p.xp()  = sin(u);
		u = RandomNG::uniform(-pi, pi);
		p.y()   = cos(u);
		p.yp()  = sin(u);
		p.dp()  = RandomNG::uniform(-1, 1);
		p.ct()  = RandomNG::uniform(-1, 1);
		break;
	case horizontalHaloDistribution1:
		u = RandomNG::uniform(-pi, pi);
		p.x()   = cos(u);
		p.xp()  = sin(u);
		u = RandomNG::uniform(-pi, pi);
		p.y()   = 0.0;
		p.yp()  = 0.0;
		p.dp()  = RandomNG::uniform(-1, 1);
		p.ct()  = RandomNG::uniform(-1, 1);
		break;
	case verticalHaloDistribution1:
		u = RandomNG::uniform(-pi, pi);
		p.x()   = 0.0;
		p.xp()  = 0.0;
		u = RandomNG::uniform(-pi, pi);
		p.y()   = cos(u);
		p.yp()  = sin(u);
		p.dp()  = RandomNG::uniform(-1, 1);
		p.ct()  = RandomNG::uniform(-1, 1);
		break;
	case horizontalHaloDistribution2:
		u = RandomNG::uniform(-pi, pi);
		p.x()   = cos(u);
		p.xp()  = sin(u);
		u = RandomNG::uniform(-pi, pi);
		p.y()   = RandomGauss(1, cutoffs.y());
		p.yp()  = RandomGauss(1, cutoffs.yp());
		p.dp()  = RandomNG::uniform(-1, 1);
		p.ct()  = RandomNG::uniform(-1, 1);
		break;
	case verticalHaloDistribution2:
		u = RandomNG::uniform(-pi, pi);
		p.x()   = RandomGauss(1, cutoffs.x());
		p.xp()  = RandomGauss(1, cutoffs.xp());
		u = RandomNG::uniform(-pi, pi);
		p.y()   = cos(u);
		p.yp()  = sin(u);
		p.dp()  = RandomNG::uniform(-1, 1);
		p.ct()  = RandomNG::uniform(-1, 1);
		break;
	}
	return p;
}

Bunch* ParticleBunchConstructor::ConstructBunch(int bunchIndex) const
{
	ParticleBunchConstructor::ConstructBunchDistribution();
	return new ParticleBunch(beamdat.p0, beamdat.charge, pbunch);
}

void ParticleBunchConstructor::ForceCentroid(bool fc)
{
	force_c = fc;
}

void ParticleBunchConstructor::DoForceCentroid() const
{
	PSvector xm = pbunch.front();
	for(auto p = pbunch.begin() + 1; p != pbunch.end(); ++p)
	{
		xm += *p;
	}

	xm /= np;
	xm -= pbunch.front();

	for(auto p = pbunch.begin() + 1; p != pbunch.end(); ++p)
	{
		p->x() -= xm.x();
		p->xp() -= xm.xp();
		p->y() -= xm.y();
		p->yp() -= xm.yp();
		p->dp() -= xm.dp();
		p->ct() -= xm.ct();
	}
}

} //end namespace ParticleTracking
