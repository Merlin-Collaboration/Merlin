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

#ifndef BetatronTunes_h
#define BetatronTunes_h 1

#include "AcceleratorModel/AcceleratorModel.h"

#include "BeamDynamics/ParticleTracking/HollowELensProcess.h"

#include "BeamModel/PSTypes.h"

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

	HollowELensProcess* myHELProcess;
	void SetHELProcess(HollowELensProcess* HELP)
	{
		myHELProcess = HELP;
	}

	double FindTune(vector<double>& data);

private:
	AcceleratorModel* theModel;
	double p0;
	void FFT(vector<double>& data);
	double amp(double a, double b, double c);
};

#endif
