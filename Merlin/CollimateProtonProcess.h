/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CollimateProtonProcess_h
#define CollimateProtonProcess_h 1

#include <fstream>

#include "CollimateParticleProcess.h"
#include "ScatteringModel.h"

using namespace Collimation;

namespace ParticleTracking
{

class CollimateProtonProcess: public CollimateParticleProcess
{
public:

	/**
	 * Constructor taking the collimation mode, and the output
	 * stream pointer to which to print the results. mode can
	 * be a logical OR combination of the collimation modes. A
	 * null pointer for osp (default) suppresses output.
	 */
	CollimateProtonProcess(int priority, int mode, std::ostream* osp = nullptr);

	void SetScatteringModel(Collimation::ScatteringModel* s);

private:
	Collimation::ScatteringModel* scattermodel;

	bool DoScatter(Particle&);

};

} // end namespace ParticleTracking

#endif
