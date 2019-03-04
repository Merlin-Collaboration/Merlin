/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ConstantGradientAdjustment
#define _h_ConstantGradientAdjustment

#include <vector>
#include "EnergyAdjustmentPolicy.h"

// Produces an arbitrary number of energy states by uniformly adjusting
// the entire gradient of the machine.
class ConstantGradientAdjustment: public EnergyAdjustmentPolicy
{
public:

	ConstantGradientAdjustment();
	~ConstantGradientAdjustment();

	size_t GetNumEnergyStates()
	{
		return estates.size();
	}

	void Initialise();
	void SetEnergyState(size_t nes);
	void SetActiveBeamlineSegment(DFS_Segment &seg)
	{
		// noting to do;
	}
	bool SupportsIncrementalTracking()
	{
		return true;
	}

	// Energy state definition:
	// Add an energy state defined by a relative
	// uniform change in accelerator gradient
	// (kystron voltage, dEk) and a relative change
	// in energy of the initial beam (dEb). Energy
	// states are added in the order that they will
	// be applied in the algorithm.
	void AddEnergyState(double dEk, double dEb);

private:

	struct EnergyState
	{
		double dEkly;
		double dEbeam;
		EnergyState(double dEk, double dEb) :
			dEkly(dEk), dEbeam(dEb)
		{
		}

	};

	std::vector<EnergyState> estates;
};

#endif
