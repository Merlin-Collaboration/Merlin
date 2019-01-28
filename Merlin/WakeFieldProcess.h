/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_WakeFieldProcess
#define _h_WakeFieldProcess

#ifdef ENABLE_MPI
#include <time.h>
#endif

#include <vector>
#include <typeinfo>

#include "WakePotentials.h"
#include "ParticleBunchProcess.h"
#include "StringPattern.h"

class WakePotentials;

namespace ParticleTracking
{

/**
 * class WakefieldProcess
 * Implements single-bunch wakefields.
 *
 * Modified by A.Wolski 12/2/2003
 * to include the derivative of the charge distribution
 * for handling CSR wake fields.
 *
 * Modified by D.Kruecker 18.2.2008
 * to be used as base class for other wakefield types (collimator,coupler,...)
 */
class WakeFieldProcess: public ParticleBunchProcess
{
public:

	enum ImpulseLocation
	{
		atCentre,
		atExit

	};

	WakeFieldProcess(int prio, size_t nb = 100, double ns = 3.0, string aID = "WAKEFIELD");

	~WakeFieldProcess();

	virtual void InitialiseProcess(Bunch& bunch);
	virtual void SetCurrentComponent(AcceleratorComponent& component);
	virtual void DoProcess(double ds);
	virtual double GetMaxAllowedStepSize() const;

	void ApplyImpulseAt(ImpulseLocation loc)
	{
		imploc = loc;
	}
	void IncludeTransverseWake(bool flg)
	{
		inc_tw = flg;
	}
	void DumpSliceCentroids(ostream&) const;
	void SetFilter(int n, int m, int d);

protected:

	ImpulseLocation imploc;
	double current_s;
	double impulse_s;
	double clen;
	size_t nbins;
	double nsig;

	void Init();
	size_t CalculateQdist();
	virtual void CalculateWakeL();
	virtual void CalculateWakeT();
	virtual void ApplyWakefield(double ds);

	WakePotentials* currentWake;

	std::vector<ParticleBunch::iterator> bunchSlices;
	std::vector<double> Qd;
	std::vector<double> Qdp;
	std::vector<double>* filter;
	std::vector<double> wake_x;
	std::vector<double> wake_y;
	std::vector<double> wake_z;
	bool recalc;
	bool inc_tw;

	double zmin, zmax, dz;

	size_t oldBunchLen;

private:
	//Copy protection
	WakeFieldProcess(const WakeFieldProcess& rhs);
	WakeFieldProcess& operator=(const WakeFieldProcess& rhs);
};

void savgol(vector<double>& c, int nl, int nr, int ld, int m);

} // end namespace ParticleTracking

#endif
