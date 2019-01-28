/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ParticleBunch.h"
#include "ParticleTracker.h"
#include "SynchRadParticleProcess.h"
#include "RingDeltaTProcess.h"
#include "ClosedOrbit.h"
#include "TransferMatrix.h"
#include "MatrixPrinter.h"

using namespace ParticleTracking;

TransferMatrix::TransferMatrix(AcceleratorModel* aModel, double refMomentum) :
	theModel(aModel), p0(refMomentum), radiation(false), obspnt(0), delta(1.0e-9), bendscale(0)
{
}

void TransferMatrix::Radiation(bool flag)
{
	radiation = flag;

	if(radiation)
	{
		SetRadNumSteps(1);
	}
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
	{
		co.SetRadNumSteps(radnumsteps);
	}
	else
	{
		co.SetRadStepSize(radstepsize);
	}

	if(bendscale != 0)
	{
		co.ScaleBendPathLength(bendscale);
	}

	co.FindClosedOrbit(p, obspnt);
	FindTM(M, p);
}

void TransferMatrix::FindClosedOrbitTM(RealMatrix& M, PSvector& orbit)
{
	ClosedOrbit co(theModel, p0);
	co.Radiation(radiation);

	if(radstepsize == 0)
	{
		co.SetRadNumSteps(radnumsteps);
	}
	else
	{
		co.SetRadStepSize(radstepsize);
	}

	if(bendscale != 0)
	{
		co.ScaleBendPathLength(bendscale);
	}

	co.FindClosedOrbit(orbit, obspnt);
	FindTM(M, orbit);
}

void TransferMatrix::FindTM(RealMatrix& M, PSvector& orbit)
{
	ParticleBunch bunch(p0, 1.0);
	int k = 0;
	for(k = 0; k < 7; k++)
	{
		Particle p = orbit;
		if(k > 0)
		{
			p[k - 1] += delta;
		}
		bunch.push_back(p);
	}
//cout <<"bunch dump 1" << endl;
	//MatrixForm(M,std::cout,OPFormat().precision(6).fixed());
//	bunch.Output(std::cout);
	ParticleTracker tracker(theModel->GetRing(obspnt), &bunch, false);

	if(radiation)
	{
		SynchRadParticleProcess* srproc = new SynchRadParticleProcess(1);

		if(radstepsize == 0)
		{
			srproc->SetNumComponentSteps(radnumsteps);
		}
		else
		{
			srproc->SetMaxComponentStepSize(radstepsize);
		}

		srproc->AdjustBunchReferenceEnergy(false);
		tracker.AddProcess(srproc);
	}

	if(bendscale != 0)
	{
		RingDeltaTProcess* ringdt = new RingDeltaTProcess(2);
		ringdt->SetBendScale(bendscale);
		tracker.AddProcess(ringdt);
	}

	tracker.Run();

	ParticleBunch::const_iterator ip = tracker.GetTrackedBunch().begin();
	const Particle& pref = *ip++;

	for(k = 0; k < 6; k++, ip++)
		for(int m = 0; m < 6; m++)
		{
			M(m, k) = ((*ip)[m] - pref[m]) / delta;
		}

//      cout<<"Transfer Matrix: "<<endl;
//      cout<<orbit;
//      cout<<pref<<endl;
//      MatrixForm(M,cout,OPFormat().precision(6).fixed());
//      cout<<endl;
}

void TransferMatrix::FindTM(RealMatrix& M, PSvector& orbit, int n1, int n2)
{
	ParticleBunch bunch(p0, 1.0);
	int k = 0;
	for(k = 0; k < 7; k++)
	{
		Particle p = orbit;
		if(k > 0)
		{
			p[k - 1] += delta;
		}
		bunch.push_back(p);
	}

	ParticleTracker tracker(theModel->GetBeamline(n1, n2), &bunch, false);

	if(radiation)
	{
		SynchRadParticleProcess* srproc = new SynchRadParticleProcess(1);

		if(radstepsize == 0)
		{
			srproc->SetNumComponentSteps(radnumsteps);
		}
		else
		{
			srproc->SetMaxComponentStepSize(radstepsize);
		}

		srproc->AdjustBunchReferenceEnergy(false);
		tracker.AddProcess(srproc);
	}

	if(bendscale != 0)
	{
		RingDeltaTProcess* ringdt = new RingDeltaTProcess(2);
		ringdt->SetBendScale(bendscale);
		tracker.AddProcess(ringdt);
	}

	tracker.Run();

	ParticleBunch::const_iterator ip = tracker.GetTrackedBunch().begin();
	const Particle& pref = *ip++;

	for(k = 0; k < 6; k++, ip++)
		for(int m = 0; m < 6; m++)
		{
			M(m, k) = ((*ip)[m] - pref[m]) / delta;
		}
}
