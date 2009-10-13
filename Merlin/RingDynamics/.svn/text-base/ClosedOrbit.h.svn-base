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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ClosedOrbit_h
#define ClosedOrbit_h 1

#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamDynamics/ParticleTracking/ParticleTracker.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "BeamModel/PSTypes.h"

using namespace ParticleTracking;

class ClosedOrbit
{
public:
    ClosedOrbit(AcceleratorModel* aModel, double refMomentum);
    ~ClosedOrbit();

    void FindClosedOrbit(PSvector& particle, int ncpt = 0);
    void FindRMSOrbit(PSvector& particle);

    void TransverseOnly(bool flag);				// default: false
    void Radiation(bool flag);					// default: false
    void SetRadStepSize(double rad_stepsize);
    void SetRadNumSteps(int rad_numsteps);
    void FullAcceleration(bool flag);			// default: false
    void ScaleBendPathLength(double scale);

    void SetDelta(double new_delta);			// default: 1.0e-9
    void SetTolerance(double tolerance);		// default: 1.0e-26
    void SetMaxIterations(int max_iterations);	// default: 20

    void AddProcess(ParticleBunchProcess* aProcess);

    // The following member functions are available for diagnostics

    // The final achieved figure of merit for the iteration
    double w;

    // The number of iterations
    int iter;

private:
    AcceleratorModel* theModel;
    ParticleTracker* theTracker;
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
};

#endif

