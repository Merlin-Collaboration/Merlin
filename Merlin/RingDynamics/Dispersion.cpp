/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/SynchRadParticleProcess.h"
#include "RingDynamics/ClosedOrbit.h"
#include "RingDynamics/Dispersion.h"
#include "NumericalUtils/utils.h"

using namespace ParticleTracking;

Dispersion::Dispersion(AcceleratorModel* aModel, double refMomentum)
        : theModel(aModel), p0(refMomentum), delta(1.0e-9) {}

double Dispersion::SetDelta(double new_delta)
{
    double old = delta;
    delta = new_delta;
    return old;
}

void Dispersion::FindDispersion(int n)
{
    ClosedOrbit co(theModel, p0);
    co.TransverseOnly(true);

    PSvector p(0);
    p.dp() = -delta;
    co.FindClosedOrbit(p, n);

    PSvector q(0);
    q.dp() = delta;
    co.FindClosedOrbit(q, n);

    Dx  = (q.x() - p.x())/2/delta;
    Dxp = (q.xp() - p.xp())/2/delta;
    Dy  = (q.y() - p.y())/2/delta;
    Dyp = (q.yp() - p.yp())/2/delta;
}

void Dispersion::FindRMSDispersion(ofstream* file)
{
    ClosedOrbit co(theModel, p0);
    co.TransverseOnly(true);
    PSvector p(0);
    p.dp() = -delta;
    co.FindClosedOrbit(p);

    PSvector q(0);
    q.dp() = delta;
    co.FindClosedOrbit(q);

    ParticleBunch* particle = new ParticleBunch(p0, 1.0);
    particle->push_back(p);
    particle->push_back(q);

    ParticleTracker tracker(theModel->GetBeamline(),particle);

    double len = 0.0;
    double dl = 0.0;
    double prev[2];
    prev[0] = prev[1] = 0;
    double pres[2];
    double rmsD[2];
    rmsD[0] = rmsD[1] = 0;

    bool loop = true;
    tracker.InitStepper();

    do {
        if(file)
            *file<<std::setw(14)<<len;

        dl = tracker.GetCurrentComponent().GetLength();
        loop = tracker.StepComponent();

        // This is a little clumsy, but is a direct result of not
        // having random access for std::list!
        ParticleBunch::const_iterator ip = tracker.GetTrackedBunch().begin();
        const Particle& p0 = *ip++;
        const Particle& p1 = *ip;

        for(int m=0; m<2; m++)
        {
            pres[m] = (p1[2*m] - p0[2*m])/2/delta;
            rmsD[m] += dl * (pres[m]*pres[m] + prev[m]*prev[m]) / 2;
            prev[m] = pres[m];
            if(file)
                *file<<std::setw(14)<<pres[m];
        }

        if(file)
            *file<<endl;

        len += dl;

    } while(loop);

    DxRMS = sqrt(rmsD[0] / len);
    DyRMS = sqrt(rmsD[1] / len);
}
