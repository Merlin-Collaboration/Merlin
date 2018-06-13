/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_BeamDynamicsModel
#define _h_BeamDynamicsModel 1

#include "Bunch.h"
#include "BeamData.h"
#include "AcceleratorModel.h"

class SimulationOutput;

// Encapsulates the exact physics models for simulating the transport
// (acceleration) of a beam.
class BeamDynamicsModel
{
public:

	BeamDynamicsModel(const std::string& aName) :
		output(nullptr), itsName(aName)
	{
	}

	virtual ~BeamDynamicsModel()
	{
	}

	// Sets the beamline to be tracked by all subsequent
	// TrackBunch() operations.
	virtual void SetBeamline(const AcceleratorModel::Beamline& bl) = 0;

	// Sets the initial bunch to be tracked by all subsequent
	// TrackBunch() operations.
	virtual void SetInitialBunch(const Bunch* bunch0) = 0;

	// Tracks a copy of the current initial bunch through
	// the current beamline. Returns the resulting tracked
	// bunch;
	virtual Bunch* TrackBunch() = 0;

	// Tracks and updates the specified bunch. This method
	// does not affect the current initial bunch set by SetInitialBunch().
	virtual void TrackThisBunch(Bunch* b) = 0;

	// Creates a bunch with the given initial beam specification.
	virtual Bunch* CreateBunch(const BeamData& beam0) = 0;

	// Return the name of this model
	const std::string& GetName() const
	{
		return itsName;
	}

	// Set the output object to be used during tracking
	void SetOutput(SimulationOutput* so);

protected:

	SimulationOutput* output;

private:

	const std::string itsName;

};

inline void BeamDynamicsModel::SetOutput(SimulationOutput* so)
{
	output = so;
}

#endif
