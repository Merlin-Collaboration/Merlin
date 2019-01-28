/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef HollowELensProcess_h
#define HollowELensProcess_h 1

#include "merlin_config.h"
#include "AcceleratorModel.h"
#include "HollowElectronLens.h"
#include "ParticleBunchProcess.h"
#include "ParticleBunch.h"

namespace ParticleTracking
{

class HollowELensProcess: public ParticleBunchProcess
{
public:
	/**
	 *	Constructor
	 */
	HollowELensProcess(int priority);

	/**
	 *	Initialise this process with the specified Bunch. If
	 *	bunch is not a ParticleBunch object, the process becomes
	 *	inactive.
	 */
	virtual void InitialiseProcess(Bunch& bunch);

	/**
	 *	Sets the current accelerator component.
	 */
	virtual void SetCurrentComponent(AcceleratorComponent& component);

	/**
	 *	Preform the process for the specified step ds.
	 */
	virtual void DoProcess(double ds);

	/**
	 *	Returns the current maximum step length for this process.
	 *	@return Current process maximum step length
	 */
	virtual double GetMaxAllowedStepSize() const;

	/**
	 * Calculates the theta kick given by the e- lens
	 */
	virtual double CalcThetaMax(double r);

	/**
	 * Use simple profile to calculate kick
	 */
	virtual double CalcKickSimple(Particle &p);
	/**
	 * Need this to output profiles
	 */
	virtual double CalcKickSimple(double R);

	/**
	 * Use radial (measured) profile to calculate kick
	 */
	virtual double CalcKickRadial(Particle &p);
	/**
	 * Need this to output profiles
	 */
	virtual double CalcKickRadial(double R);

	/**
	 * Output the HEL radial profile in x y phase space (assumes circular HEL)
	 */
	virtual void OutputProfile(std::ostream* os, double E = 7000, double min = 0, double max = 10);

private:
	// Data Members for Class Attributes
	HollowElectronLens* currentComponentHEL;
	double ProtonBeta;
};

} // end namespace ParticleTracking
#endif
