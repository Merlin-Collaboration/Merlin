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

namespace ParticleTracking
{

class HollowELensProcess : public ParticleBunchProcess
{
public:
	/**
	*	Constructor
	*/
	HollowELensProcess (int priority);

	/**
	*	Initialise this process with the specified Bunch. If
	*	bunch is not a ParticleBunch object, the process becomes
	*	inactive.
	*/
	virtual void InitialiseProcess (Bunch& bunch);

	/**
	*	Sets the current accelerator component.
	*/
	virtual void SetCurrentComponent (AcceleratorComponent& component);

	/**
	*	Preform the process for the specified step ds.
	*/
	virtual void DoProcess (double ds);

	/**
	*	Returns the current maximum step length for this process.
	*	@return Current process maximum step length
	*/
	virtual double GetMaxAllowedStepSize () const;

	/**
	* Calculates the theta kick given by the e- lens
	*/
	virtual double CalcThetaMax (double r);

	/**
	* Use simple profile to calculate kick
	*/
	virtual double CalcKickSimple (Particle &p);
	/**
	* Need this to output profiles
	*/
	virtual double CalcKickSimple (double R);

	/**
	* Use radial (measured) profile to calculate kick
	*/
	virtual double CalcKickRadial (Particle &p);
	/**
	* Need this to output profiles
	*/
	virtual double CalcKickRadial (double R);

	/**
	* Output the HEL radial profile in x y phase space (assumes circular HEL)
	*/
	virtual void OutputProfile(std::ostream* os, double E=7000, double min=0, double max=10);

private:
	// Data Members for Class Attributes
	HollowElectronLens* currentComponentHEL;
	double ProtonBeta;
};

} // end namespace ParticleTracking
#endif
