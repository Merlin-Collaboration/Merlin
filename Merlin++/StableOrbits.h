/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef StableOrbits_h
#define StableOrbits_h 1

#include <list>
#include "AcceleratorModel.h"
#include "ParticleBunch.h"

using namespace ParticleTracking;

class StableOrbits
{
public:
	StableOrbits(AcceleratorModel* aModel);
	void SelectStable(ParticleBunch& aBunch, std::list<size_t>* index);

	int SetTurns(int turns);
	int SetObservationPoint(int n);

private:
	AcceleratorModel* theModel;
	int nturns;
	int obspnt;
};

#endif
