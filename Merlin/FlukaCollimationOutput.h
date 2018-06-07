/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef FlukaCollimationOutput_h
#define FlukaCollimationOutput_h 1

#include <string>
#include <vector>

#include "CollimationOutput.h"
#include "AcceleratorComponent.h"
#include "ParticleBunch.h"
#include "PSTypes.h"

namespace ParticleTracking
{

class FlukaCollimationOutput: public CollimationOutput
{

public:

	FlukaCollimationOutput(OutputType otype = tencm);
	~FlukaCollimationOutput();

	virtual void Finalise();
	virtual void Output(std::ostream* os);
	virtual void Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn = 0);

protected:

private:

};

}

#endif
