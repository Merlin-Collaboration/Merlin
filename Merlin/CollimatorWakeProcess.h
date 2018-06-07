/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_CollimatorWakeProcess
#define _h_CollimatorWakeProcess

#include <vector>

#include "merlin_config.h"
#include "CollimatorWakePotentials.h"
#include "ParticleBunch.h"
#include "ParticleBunchProcess.h"
#include "WakeFieldProcess.h"
#include "StringPattern.h"

namespace ParticleTracking
{

/**
 * Class for calculating the longitudinal and
 * transverse single-bunch wakefields
 * for Collimators with modes
 */
class CollimatorWakeProcess: public WakeFieldProcess
{
public:

	CollimatorWakeProcess(int, int, size_t, double);
	~CollimatorWakeProcess();

	virtual void ApplyWakefield(double);
	virtual void CalculateWakeT(double, int);
	virtual void CalculateWakeL(double, int);

private:

	double CalculateSm(int, int);
	double CalculateCm(int, int);

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

	using WakeFieldProcess::CalculateWakeT;
	using WakeFieldProcess::CalculateWakeL;
};

}

#endif
