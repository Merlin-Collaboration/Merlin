/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/04/14 16:01:43 $
// $Revision: 1.9 $
//
/////////////////////////////////////////////////////////////////////////
#include <memory>

#include "RandomNG.h"
#include "ParticleBunchConstructor.h"
#include "NormalTransform.h"
#include "NumericalConstants.h"

namespace ParticleTracking
{

inline double RandomGauss(double variance, double cutoff)
{
	return cutoff==0 ? RandomNG::normal(0,variance) :  RandomNG::normal(0,variance,cutoff);
}

ParticleBunchConstructor::ParticleBunchConstructor (const BeamData& beam, size_t npart, DistributionType dist)
	: np(npart),dtype(dist),cutoffs(0),beamdat(beam),itsFilter(nullptr),M(NormalTransform(beam)),force_c(false)
{}

ParticleBunchConstructor::~ParticleBunchConstructor ()
{
	if(itsFilter)
	{
		delete itsFilter;
	}
}

void ParticleBunchConstructor::SetBunchData (const BeamData& beam)
{
	beamdat = beam;
	M.R = NormalTransform(beam);
}

void ParticleBunchConstructor::SetNumParticles (size_t npart)
{
	assert(npart>0);
	np=npart;
}

void ParticleBunchConstructor::SetDistributionCutoff (double cut)
{
	cutoffs = PSvector(fabs(cut));
}

void ParticleBunchConstructor::SetDistributionCutoff (const PSvector& cut)
{
	cutoffs=cut;
}

void ParticleBunchConstructor::ConstructBunchDistribution (int bunchIndex) const
{


	// The first particle is *always* the centroid particle
	PSvector p;
	p.x()=beamdat.x0;
	p.xp()=beamdat.xp0;
	p.y()=beamdat.y0;
	p.yp()=beamdat.yp0;
	p.dp()=0;
	p.ct()=beamdat.ct0;
	p.type() = -1.0;
	p.location() = -1.0;
	p.id() = 0;
	p.sd() = 0.0;
	pbunch.push_back(p);

	size_t i = 1;
	while(i<np)
	{
		p = GenerateFromDistribution();
		M.Apply(p);
		p+=pbunch.front(); // add centroid
		p.type() = -1.0;
		p.location() = -1.0;
		p.id() = i;
		p.sd() = 0.0;
		if(itsFilter==nullptr || itsFilter->Apply(p))
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
	double rx, ry, u;
	switch(dtype)
	{
	case normalDistribution:
		p.x()	= RandomGauss(beamdat.emit_x,cutoffs.x());
		p.xp()	= RandomGauss(beamdat.emit_x,cutoffs.xp());
		p.y()	= RandomGauss(beamdat.emit_y,cutoffs.y());
		p.yp()	= RandomGauss(beamdat.emit_y,cutoffs.yp());
		p.dp()	= RandomGauss(beamdat.sig_dp*beamdat.sig_dp,cutoffs.dp());
		p.ct()	= RandomGauss(beamdat.sig_z*beamdat.sig_z,cutoffs.ct());
		break;
	case flatDistribution:
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		p.x()	= RandomNG::uniform(-rx,rx);
		p.xp()	= RandomNG::uniform(-rx,rx);
		p.y()	= RandomNG::uniform(-ry,ry);
		p.yp()	= RandomNG::uniform(-ry,ry);
		p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
		p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
		break;
	case skewHaloDistribution:
	case ringDistribution:
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		u = RandomNG::uniform(-pi,pi);
		p.x()	= rx * cos(u);
		p.xp()	= rx * sin(u);
		u = RandomNG::uniform(-pi,pi);
		p.y()	= ry * cos(u);
		p.yp()	= ry * sin(u);
		p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
		p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
		break;
	case horizontalHaloDistribution1:
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		u = RandomNG::uniform(-pi,pi);
		p.x()	= rx * cos(u);
		p.xp()	= rx * sin(u);
		u = RandomNG::uniform(-pi,pi);
		p.y()	= 0.0;
		p.yp()	= 0.0;
		p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
		p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
		break;
	case verticalHaloDistribution1:
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		u = RandomNG::uniform(-pi,pi);
		p.x()	= 0.0;
		p.xp()	= 0.0;
		u = RandomNG::uniform(-pi,pi);
		p.y()	= ry * cos(u);
		p.yp()	= ry * sin(u);
		p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
		p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
		break;
	case horizontalHaloDistribution2:
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		u = RandomNG::uniform(-pi,pi);
		p.x()	= rx * cos(u);
		p.xp()	= rx * sin(u);
		u = RandomNG::uniform(-pi,pi);
		p.y()	= RandomGauss(beamdat.emit_y,cutoffs.y());
		p.yp()	= RandomGauss(beamdat.emit_y,cutoffs.yp());
		p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
		p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
		break;
	case verticalHaloDistribution2:
		rx = sqrt(beamdat.emit_x);
		ry = sqrt(beamdat.emit_y);
		u = RandomNG::uniform(-pi,pi);
		p.x()	= RandomGauss(beamdat.emit_x,cutoffs.x());
		p.xp()	= RandomGauss(beamdat.emit_x,cutoffs.xp());
		u = RandomNG::uniform(-pi,pi);
		p.y()	= ry * cos(u);
		p.yp()	= ry * sin(u);
		p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
		p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
		break;
	}
	return p;
}

Bunch* ParticleBunchConstructor::ConstructBunch (int bunchIndex) const
{
	ParticleBunchConstructor::ConstructBunchDistribution();
	return new ParticleBunch(beamdat.p0,beamdat.charge,pbunch);
}

void ParticleBunchConstructor::ForceCentroid (bool fc)
{
	force_c = fc;
}

void ParticleBunchConstructor::DoForceCentroid () const
{
	PSvector xm = pbunch.front();
	for (auto p = pbunch.begin()+1; p != pbunch.end(); ++p)
	{
		xm += *p;
	}

	xm /= np;
	xm -= pbunch.front();

	for (auto p = pbunch.begin()+1; p != pbunch.end(); ++p)
	{
		(*p) -= xm;
	}
}


} //end namespace ParticleTracking
