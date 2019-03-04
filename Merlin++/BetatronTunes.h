/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef BetatronTunes_h
#define BetatronTunes_h 1

#include "AcceleratorModel.h"
#include "HollowELensProcess.h"
#include "PSTypes.h"

using namespace ParticleTracking;

class BetatronTunes
{
public:
	BetatronTunes(AcceleratorModel* aModel, double refMomentum);
	void FindTunes(PSvector& particle, int ntrack = 256, bool diffusion = true);
	double Qx, Qy, dQx, dQy;

	double GetQx()
	{
		return Qx;
	}
	double GetdQx()
	{
		return dQx;
	}
	double GetQy()
	{
		return Qy;
	}
	double GetdQy()
	{
		return dQy;
	}

	void SetHELProcess(HollowELensProcess* HELP)
	{
		myHELProcess = HELP;
	}

	double FindTune(vector<double>& data);

private:
	AcceleratorModel* theModel;
	double p0;
	HollowELensProcess* myHELProcess;

	void FFT(vector<double>& data);
	double amp(double a, double b, double c);
};

#endif
