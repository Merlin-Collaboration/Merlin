/////////////////////////////////////////////////////////////////////////
// class KSEnergyAdjustment
// Klystron Shunting Energy Adjustment performs a LIAR-like energy
// modification by turning off a set of klystrons upstream of the
// current correction segment. The number of klystrons to turn off
// is calculated based on the required relative energy difference.
// If not enough klystron are available upstream of the segment (as is
// the case for the first segment for example), then the remaining
// energy difference is made up by adjusting the initial beam energy.
//
// Note that this method does not support incremental tracking.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/19 10:19:06 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_KSEnergyAdjustment
#define _h_KSEnergyAdjustment 1

#include "EnergyAdjustmentPolicy.h"

class KSEnergyAdjustment :public EnergyAdjustmentPolicy {

public:

	explicit KSEnergyAdjustment(double dEr);

	// Only the nominal beam and a single off-energy
	// state are current supported.
	size_t GetNumEnergyStates() {
		return 2; 
	}

	// Sets the nes-th energy state.
	void SetEnergyState(size_t nes);

	// Initialise the policy. This function is called 
	// before implementation of DFS to allow any implementation
	// dependent initialisation to be performed.
	void Initialise();

	// Informs this policy of the beamline segment being tuned.
	void SetActiveBeamlineSegment(DFS_Segment &seg);

	// Incremental tracking is not supported.
	bool SupportsIncrementalTracking() {
		return false;
	}

private:

	double energy0; 
	double delta;
	std::pair<size_t,size_t> cKlysRange;
	double dEbeam;
};

#endif

