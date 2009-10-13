#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/SynchRadParticleProcess.h"
#include "BeamDynamics/ParticleTracking/RingDeltaTProcess.h"
#include "RingDynamics/ClosedOrbit.h"
#include "RingDynamics/TransferMatrix.h"
#include "NumericalUtils/MatrixPrinter.h"

using namespace ParticleTracking;

TransferMatrix::TransferMatrix(AcceleratorModel* aModel, double refMomentum)
        : theModel(aModel), p0(refMomentum), radiation(false), obspnt(0), delta(1.0e-9), bendscale(0) {}

void TransferMatrix::Radiation(bool flag)
{
    radiation = flag;

    if(radiation)
        SetRadNumSteps(1);
}

void TransferMatrix::SetObservationPoint(int n)
{
    obspnt = n;
}
void TransferMatrix::SetDelta(double new_delta)
{
    delta = new_delta;
}

void TransferMatrix::SetRadStepSize(double rad_stepsize)
{
    radstepsize = rad_stepsize;
    radnumsteps = 0;
}

void TransferMatrix::SetRadNumSteps(int rad_numsteps)
{
    radnumsteps = rad_numsteps;
    radstepsize = 0;
}

void TransferMatrix::ScaleBendPathLength(double scale)
{
    bendscale = scale;
}

void TransferMatrix::FindTM(RealMatrix& M)
{
    PSvector p(0);
    ClosedOrbit co(theModel, p0);
    co.Radiation(radiation);

    if(radstepsize == 0)
        co.SetRadNumSteps(radnumsteps);
    else
        co.SetRadStepSize(radstepsize);

    if(bendscale!=0)
        co.ScaleBendPathLength(bendscale);

    co.FindClosedOrbit(p, obspnt);
    FindTM(M,p);
}

void TransferMatrix::FindClosedOrbitTM(RealMatrix& M, PSvector& orbit)
{
    ClosedOrbit co(theModel, p0);
    co.Radiation(radiation);

    if(radstepsize == 0)
        co.SetRadNumSteps(radnumsteps);
    else
        co.SetRadStepSize(radstepsize);

    if(bendscale!=0)
        co.ScaleBendPathLength(bendscale);

    co.FindClosedOrbit(orbit, obspnt);
    FindTM(M,orbit);
}

void TransferMatrix::FindTM(RealMatrix& M, PSvector& orbit)
{
    ParticleBunch bunch(p0,1.0);
    int k=0;
    for(k=0; k<7; k++)	{
        Particle p = orbit;
        if(k>0)
            p[k-1] += delta;
        bunch.push_back(p);
    }

    ParticleTracker tracker(theModel->GetRing(obspnt),&bunch,false);

    if(radiation) {
        SynchRadParticleProcess* srproc = new SynchRadParticleProcess(1);

        if(radstepsize == 0)
            srproc->SetNumComponentSteps(radnumsteps);
        else
            srproc->SetMaxComponentStepSize(radstepsize);

        srproc->AdjustBunchReferenceEnergy(false);
        tracker.AddProcess(srproc);
    }

    if(bendscale!=0) {
        RingDeltaTProcess* ringdt = new RingDeltaTProcess(2);
        ringdt->SetBendScale(bendscale);
        tracker.AddProcess(ringdt);
    }

    tracker.Run();

    ParticleBunch::const_iterator ip = tracker.GetTrackedBunch().begin();
    const Particle& pref = *ip++;

    for(k=0; k<6; k++,ip++)
        for(int m=0; m<6; m++)
            M(m,k) = ((*ip)[m] - pref[m]) / delta;

    /*
    	cout<<"Transfer Matrix: "<<endl;
    	cout<<orbit;
    	cout<<pref<<endl;
    	MatrixForm(M,cout,OPFormat().precision(6).fixed());
    	cout<<endl;
    */
}
