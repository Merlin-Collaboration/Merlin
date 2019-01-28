/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <fstream>

#include "ParticleTracker.h"
#include "CollimateParticleProcess.h"
#include "StableOrbits.h"

StableOrbits::StableOrbits(AcceleratorModel* aModel) :
	theModel(aModel), nturns(1), obspnt(0)
{
}

int StableOrbits::SetTurns(int turns)
{
	int old = nturns;
	nturns = turns;
	return old;
}

int StableOrbits::SetObservationPoint(int n)
{
	int old = obspnt;
	obspnt = n;
	return old;
}

void StableOrbits::SelectStable(ParticleBunch& bunch, list<size_t>* index)
{
	ParticleTracker tracker(theModel->GetRing(obspnt), &bunch, false);

	ofstream collimlog("DataFiles/CollimationLog.dat");
	CollimateParticleProcess* collimate = new CollimateParticleProcess(1, COLL_AT_CENTER, &collimlog);
	collimate->IndexParticles(*index);
	tracker.AddProcess(collimate);

	for(int turn_count = 1; turn_count <= nturns; turn_count++)
	{
		if(turn_count == 1)
		{
			tracker.Run();
		}
		else
		{
			tracker.Continue();
		}

		ParticleBunch& tracked_bunch = tracker.GetTrackedBunch();
		cout << "Tracked turn " << turn_count << ": ";
		cout << tracked_bunch.size() << " particles remaining.    " << char(13);
	}

	cout << endl;
}
