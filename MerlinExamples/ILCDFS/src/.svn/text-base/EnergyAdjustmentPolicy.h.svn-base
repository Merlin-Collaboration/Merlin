/////////////////////////////////////////////////////////////////////////
// Abstract class EnergyAdjustmentPolicy
// Encapsulates the method for adjusting the energy of the accelerator for DFS. 
// Each energy measurement is represented by an energy state, of which
// there can be any number. Each energy state is indexed by an integer ranging
// from 1... number of energy states supported. Energy state 0 is the nominal
// design configuration of the accelerator.
//
// Energy adjustment is possible via adjusting the klystrons (physically 
// realistic method), or by adjusting the tracked beam energy directly
// (physically unrealistic).
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

#ifndef _h_EnergyAdjustmentPolicy
#define _h_EnergyAdjustmentPolicy

#include "CommonDataStructures.h"

class EnergyAdjustmentPolicy {

public:

	virtual ~EnergyAdjustmentPolicy() {}

	// Return the total number of energy states supported.
	virtual size_t GetNumEnergyStates() = 0;

	// Sets the nes-th energy state.
	virtual void SetEnergyState(size_t nes) = 0;

	// Restore the energy state to nominal value.
	virtual void Restore() {
		SetEnergyState(0);
	}

	// Initialise the policy. This function is called 
	// before implementation of DFS to allow any implementation
	// dependent initialisation to be performed.
	virtual void Initialise() =0;

	// Informs this policy of the beamline segment being tuned.
	virtual void SetActiveBeamlineSegment(DFS_Segment &seg) = 0;

	// Returns true if this policy supports incremental tracking.
	virtual bool SupportsIncrementalTracking() = 0;	

	// Attach the list of klystrons to be used for subsequent
	// adjustments.
	void SetKlystrons(const KlystronArray& klys);

	// Attach the list of ReferenceParticles (bunches) to be used for
	// subsequent adjustments.
	void SetReferenceParticles(const ReferenceParticleArray& refplist);

protected:

	KlystronArray theKlystrons;
	std::vector<Complex> defkvals;

	// Pointers to the cached bunches being tracked.
	ReferenceParticleArray beamrefs;

	// utility functions
	void RestoreKlystrons();
};

inline 	void EnergyAdjustmentPolicy::SetReferenceParticles(const ReferenceParticleArray& refplist)
{
	beamrefs = refplist;
}

#endif

