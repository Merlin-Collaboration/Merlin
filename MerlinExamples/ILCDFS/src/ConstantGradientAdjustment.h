/////////////////////////////////////////////////////////////////////////
// class ConstantGradientAdjustment
// Produces an arbitrary number of energy states by uniformly adjusting
// the entire gradient of the machine.
//
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_ConstantGradientAdjustment
#define _h_ConstantGradientAdjustment

#include <vector>
#include "EnergyAdjustmentPolicy.h"

class ConstantGradientAdjustment : public EnergyAdjustmentPolicy {	
public:

	ConstantGradientAdjustment();	
	~ConstantGradientAdjustment();	

	size_t GetNumEnergyStates() {
		return estates.size();
	}

	void Initialise();
	void SetEnergyState(size_t nes);	
	void SetActiveBeamlineSegment(DFS_Segment &seg) {
		// noting to do;
	}
	bool SupportsIncrementalTracking() { 
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

	struct EnergyState {
		double dEkly;
		double dEbeam;
		EnergyState(double dEk, double dEb) 
			: dEkly(dEk),dEbeam(dEb) {}
	};

	std::vector<EnergyState> estates;
};

#endif

