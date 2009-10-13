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

#include "BasicTransport/NormalTransform.h"
#include <memory>
// RandomNG
#include "Random/RandomNG.h"
// ParticleBunchConstructor
#include "BeamDynamics/ParticleTracking/ParticleBunchConstructor.h"

namespace ParticleTracking {

inline double RandomGauss(double variance, double cutoff)
{
    return cutoff==0 ? RandomNG::normal(0,variance) :  RandomNG::normal(0,variance,cutoff);
}

ParticleBunchConstructor::ParticleBunchConstructor (const BeamData& beam, size_t npart, DistributionType dist)

        : np(npart),dtype(dist),cutoffs(0),beamdat(beam),itsFilter(0),M(NormalTransform(beam)),force_c(false)
{}

ParticleBunchConstructor::~ParticleBunchConstructor ()
{
    if(itsFilter)
        delete itsFilter;
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

Bunch* ParticleBunchConstructor::ConstructBunch (int bunchIndex) const
{
    PSvectorArray pbunch;
    PSvector p;

    // First we generate npart particles in "normalised" phase
    // space, after which we transform them to "real" phase
    // space using M.

    // The first particle is *always* the centroid particle
    double dp2 = pow(beamdat.sig_dp,2);
    double dz2 = pow(beamdat.sig_z,2);
    double rx,ry;

    p.x()=beamdat.x0;
    p.xp()=beamdat.xp0;
    p.y()=beamdat.y0;
    p.yp()=beamdat.yp0;
    p.dp()=0;
    p.ct()=beamdat.ct0;
    pbunch.push_back(p);

    size_t i;

    PSvector xm = p; // used for calculating mean

    switch(dtype) {
    case normalDistribution:
        for(i=1; i<np;) {
            p.x()	= RandomGauss(beamdat.emit_x,cutoffs.x());
            p.xp()	= RandomGauss(beamdat.emit_x,cutoffs.xp());
            p.y()	= RandomGauss(beamdat.emit_y,cutoffs.y());
            p.yp()	= RandomGauss(beamdat.emit_y,cutoffs.yp());
            p.dp()	= RandomGauss(dp2,cutoffs.dp());
            p.ct()	= RandomGauss(dz2,cutoffs.ct());

            M.Apply(p);
            p+=pbunch.front(); // add centroid

            if(itsFilter==0 || itsFilter->Apply(p)) {
                pbunch.push_back(p);
                xm += p;
                i++;
            }

        }
        if(force_c) {
            xm/=np;
            xm-=pbunch.front();
            PSvectorArray::iterator pp=pbunch.begin();
            pp++;
            for(;pp!=pbunch.end(); pp++)
                (*pp)-=xm;
        }
        break;
    case flatDistribution:
        rx = sqrt(beamdat.emit_x);
        ry = sqrt(beamdat.emit_y);
        for(i=1; i<np;) {
            p.x()	= RandomNG::uniform(-rx,rx);
            p.xp()	= RandomNG::uniform(-rx,rx);
            p.y()	= RandomNG::uniform(-ry,ry);
            p.yp()	= RandomNG::uniform(-ry,ry);
            p.dp()	= RandomNG::uniform(-beamdat.sig_dp,beamdat.sig_dp);
            p.ct()	= RandomNG::uniform(-beamdat.sig_z,beamdat.sig_z);
            M.Apply(p);
            p+=pbunch.front(); // add centroid

            if(itsFilter==0 || itsFilter->Apply(p)) {
                pbunch.push_back(p);
                i++;
            }
        }
        break;
    };

    return new ParticleBunch(beamdat.p0,beamdat.charge,pbunch);
}

void ParticleBunchConstructor::ForceCentroid (bool fc)
{
    force_c = fc;
}

ParticleBunchFilter::~ParticleBunchFilter ()
{
    // Nothing to do
}

}; //end namespace ParticleTracking


