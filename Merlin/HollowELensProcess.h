/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 5.01 (2015)
//
// Copyright: see Merlin/copyright.txt
//
// Created:		2014 HR
// Modified:	13.03.16 HR
// Last Edited: 13.03.16 HR
//
/////////////////////////////////////////////////////////////////////////
#ifndef HollowELensProcess_h
#define HollowELensProcess_h 1

#include "merlin_config.h"

#include "AcceleratorModel.h"
#include "HollowElectronLens.h"

#include "ParticleBunchProcess.h"
#include "ParticleBunch.h"

#include "LatticeFunctions.h"

namespace ParticleTracking
{

// HEL operation modes
typedef enum {DC, AC, Diffusive, Turnskip} OperationMode;


class HollowELensProcess : public ParticleBunchProcess
{
public:
	//	Constructor
	HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e);
	HollowELensProcess (int priority, int mode, double current, double beta_e, double rigidity, double length_e, double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss);
	//	Initialise this process with the specified Bunch. If
	//	bunch is not a ParticleBunch object, the process becomes
	//	inactive.
	virtual void InitialiseProcess (Bunch& bunch);

	//	Sets the current accelerator component.
	virtual void SetCurrentComponent (AcceleratorComponent& component);

	//	Preform the process for the specified step ds.
	virtual void DoProcess (double ds);

	//	Returns the current maximum step length for this process.
	virtual double GetMaxAllowedStepSize () const;

	// Set minimum and maximum e- beam radii in [m] or [sigma]
	virtual void SetRadii (double rmin, double rmax);
	virtual void SetRadiiSigma (double rmin, double rmax, AcceleratorModel* model, double emittance_x, double emittance_y, LatticeFunctionTable* twiss, double P0=0);

	// Set the effective length of the e- lens
	virtual void SetEffectiveLength (double l_e)
	{
		EffectiveLength = l_e;
	}

	// Calculates the theta kick given by the e- lens
	virtual double CalcThetaMax (double r);

	// Set the type of HEL operation required
	virtual void SetOpMode (OperationMode mode)
	{
		OMode = mode;
	}

	// Set variables for AC mode operation
	virtual void SetAC (double tune, double deltatune, double tunevarperstep, double turnsperstep, double multi);

	//SetTurnSkip
	virtual void SetTurnskip (int skip);

	// Use simple profile to calculate kick
	virtual double CalcKickSimple (Particle &p);
	// Need this to output profiles
	virtual double CalcKickSimple (double R);

	// Use radial (measured) profile to calculate kick
	virtual double CalcKickRadial (Particle &p);
	// Need this to output profiles
	virtual double CalcKickRadial (double R);

	// Change to radial (measured) profile, simple (perfect) is default
	virtual void SetRadialProfile()
	{
		SimpleProfile = 0;
	}
	virtual void SetPerfectProfile()
	{
		SimpleProfile = 1;
	}
	virtual void SetLHCRadialProfile()
	{
		LHC_Radial = 1;
	}

	// Change electron direction (defualt opposite protons = 1)
	virtual void SetElectronDirection(bool dir);

	// Output the HEL radial profile in x y phase space (assumes circular HEL)
	virtual void OutputProfile(std::ostream* os, double E=7000, double min=0, double max=10);

	// Output the HEL footprint in x y phase space using a mapping of particles
	virtual void OutputFootprint(std::ostream* os, int npart = 1E3);

private:
	// Data Members for Class Attributes

	// Hardware parameters
	double Current;
	double ElectronBeta;
	double Rigidity;
	double ProtonBeta;
	double EffectiveLength;
	double Rmin;
	double Rmax;
	double Sigma_x;
	double Sigma_y;

	// Variables
	double XOffset;
	double YOffset;

	//For AC mode
	double Tune;
	double DeltaTune;
	double TuneVarPerStep;
	double TurnsPerStep;
	double Multiplier;
	double Nstep;

	int Turn;
	int SkipTurn;

	double MinTune;
	double MaxTune;

	bool ACSet;					// AC mode variable set?
	bool SimpleProfile;			// 1 = use perfect HEL profile, 0 = use paramaterisation of measured LHC prototype cathode profile
	bool AlignedToOrbit;		// Is the HEL aligned to the closed orbit
	bool ElectronDirection;		// 1 = opposite protons (-ve kick), 0 = same as protons (+ve kick)
	bool LHC_Radial;			// 1 = use empirically scaled measured radial profile (LHC hardware), 0 = use measured radial HEL profile

	OperationMode OMode;

};


} // end namespace ParticleTracking
#endif
