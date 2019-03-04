/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_SMP_WakeFieldProcess
#define _h_SMP_WakeFieldProcess

#include "SMPBunch.h"
#include "SMPBunchProcess.h"

#include <vector>
#include <typeinfo>

class WakePotentials;

namespace SMPTracking
{

/**
 * Applies single-bunch long. and trans. wakefields to a SliceMacroParticles
 */

class WakeFieldProcess: public SMPBunchProcess
{
public:

	enum ImpulseLocation
	{
		atCentre,
		atExit

	};

	WakeFieldProcess(int prio, double slice_width = 1.0e-6, string aID = "SMP WAKEFIELD");
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

private:

	void ApplyWakefield(double ds);

	ImpulseLocation imploc;
	double current_s;
	double impulse_s;
	double clen;

	void Init();
	void PrepLWake();
	void PrepSlices();

	std::vector<double> wake_z;
	std::vector<SMPBunch::iterator> sliceBoundaries;
	std::vector<double> slice_z;
	std::vector<double> slice_q;
	double bload;
	bool recalc;
	bool inc_tw;
	const double dz; /// slice width for binning
	WakePotentials* currentWake;

};

} // end namespace SMPTracking

#endif
