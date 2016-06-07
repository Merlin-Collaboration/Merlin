/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
// Created: June, 2006
/////////////////////////////////////////////////////////////////////////

#ifndef _h_CollimatorWakeProcess
#define _h_CollimatorWakeProcess

#include <vector>

#include "merlin_config.h"

#include "Collimators/CollimatorWakePotentials.h"

#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"

#include "utility/StringPattern.h"

namespace ParticleTracking
{

/**
* Class for calculating the longitudinal and
* transverse single-bunch wakefields
* for Collimators with modes
*/
class CollimatorWakeProcess : public WakeFieldProcess
{
public:

	CollimatorWakeProcess(int, int, size_t, double);
	~CollimatorWakeProcess ();

	virtual void ApplyWakefield(double);
	virtual void CalculateWakeT(double, int);
	virtual void CalculateWakeL(double, int);

private:

	double CalculateSm (int, int);
	double CalculateCm (int, int);

	int nmodes;

//    double Cm[5][1000];
	double** Cm;
//    double Sm[5][1000];
	double** Sm;
//	double wake_s;
//	double wake_c;

	double** wake_sl;
	double** wake_cl;
	double** wake_ct;
	double** wake_st;

	CollimatorWakePotentials* collimator_wake;

	using  WakeFieldProcess::CalculateWakeT;
	using  WakeFieldProcess::CalculateWakeL;
};

}

#endif
