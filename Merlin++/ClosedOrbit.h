/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ClosedOrbit_h
#define ClosedOrbit_h 1

#include "AcceleratorModel.h"
#include "ParticleTracker.h"
#include "ParticleBunchProcess.h"
#include "PSTypes.h"

using namespace ParticleTracking;

class ClosedOrbit
{
public:
	ClosedOrbit(AcceleratorModel* aModel, double refMomentum);
	~ClosedOrbit();

	void FindClosedOrbit(PSvector& particle, int ncpt = 0);
	void FindRMSOrbit(PSvector& particle);

	void TransverseOnly(bool flag);             // default: false
	void Radiation(bool flag);                  // default: false
	void SetRadStepSize(double rad_stepsize);
	void SetRadNumSteps(int rad_numsteps);
	void FullAcceleration(bool flag);           // default: false
	void ScaleBendPathLength(double scale);

	void SetDelta(double new_delta);            // default: 1.0e-9
	void SetTolerance(double tolerance);        // default: 1.0e-26
	void SetMaxIterations(int max_iterations);  // default: 20

	void AddProcess(ParticleBunchProcess* aProcess);

	// The following member functions are available for diagnostics

	// The final achieved figure of merit for the iteration
	double w;

	// The number of iterations
	int iter;

private:
	AcceleratorModel* theModel;
	double p0;
	bool transverseOnly;
	bool radiation;
	bool useFullAcc;

	double delta;
	double tol;
	int max_iter;
	double radstepsize;
	int radnumsteps;
	double bendscale;
	ParticleTracker* theTracker;
};

#endif
