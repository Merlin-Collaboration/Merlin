/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_DFSApp
#define _h_DFSApp 1

#include <list>
#include "DFSCorrection.h"

class AcceleratorWithErrors;
class SimulationOutput;

// Root control class for applying DFS to a beamline
class DFSApp
{
public:

	DFSApp();
	~DFSApp();

	// General (common) DFS parameters
	void SetSegmentation(size_t nquad, size_t overlap);
	void SetIteration(size_t n);
	void SetCorrectionPlane(Accelerator::Plane);
	void SetCorrectionGain(double g);
	void SetFitWeights(double w_abs, double w_diff);

	// Set the accelerator models
	void SetReferenceModel(AcceleratorModel*, BeamData*);
	void SetSimulationModel(AcceleratorModel*, BeamData*);

	// Beam dynamics models
	void SetReferenceModelBeamDynamics(BeamDynamicsModel*);
	void SetSimulationModelBeamDynamicsForDFS(BeamDynamicsModel*);
	void SetSimulationModelBeamDynamicsForEmittance(BeamDynamicsModel*);

	// Energy adjustment policy
	void SetEnergyAdjustmentPolicy(EnergyAdjustmentPolicy*);

	// One-time initialisation
	void Initialise();

	// Apply DFS to current simulation model
	void Apply();

	// Return the Simulation Model
	AcceleratorWithErrors* GetSimulationModel();

	// Return the reference model
	Accelerator* GetReferenceModel();

	// Set the output stream
	void SetOutput(SimulationOutput*);

	// Incrementatl tracking
	void UseIncrementalTracking(bool for_ref, bool for_sim);

private:

	size_t nbs;
	size_t nbo;
	size_t nit;
	Accelerator::Plane pxy;
	double gain;

	std::list<DFSCorrection*> cachedDFS;
	void ClearCachedDFS();

	BeamDynamicsModel* dfsTracker;
	BeamDynamicsModel* emittanceTracker;

	// flags to indicate if incremental tracking
	// is requested
	bool ref_useIT;
	bool sim_useIT;
};

#endif
