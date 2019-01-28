/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef RingDeltaTProcess_h
#define RingDeltaTProcess_h 1

#include "merlin_config.h"

#include "ParticleBunchProcess.h"
#include "MultipoleField.h"

namespace ParticleTracking
{

class RingDeltaTProcess: public ParticleBunchProcess
{
public:
	RingDeltaTProcess(int prio);
	virtual void SetCurrentComponent(AcceleratorComponent& component);
	virtual void DoProcess(double ds);
	virtual double GetMaxAllowedStepSize() const;
	void SetBendScale(double bendscale);
protected:
private:
	double scale;
	double dL;
	double intS;
};

} // end namespace ParticleTracking
#endif
