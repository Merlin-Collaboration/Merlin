/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_DFSCorrection
#define _h_DFSCorrection 1

#include "Channels.h"
#include "LinearAlgebra.h"
#include "CommonDataStructures.h"
#include "Accelerator.h"

class EnergyAdjustmentPolicy;
class BPMDataFilter;
class SVDMatrix;

// Applies the DFS correction to a segment of accelerator.
class DFSCorrection
{
public:

	// The following objects are used by
	// all DFSCorrection objects.
	static Accelerator* theReferenceModel;
	static Accelerator* theSimulationModel;
	static EnergyAdjustmentPolicy* theEnergyAdjustmentPolicy;

	// Weights used for the DFS fit (common to all objects)
	static double w_diff;
	static double w_abs;

	// Constructor
	DFSCorrection(const DFS_Segment& aSegment, Accelerator::Plane xy);

	// Record data and apply the correction
	// to the simulation model
	void Initialise();
	void RecordTrajectories();
	void CalculateCorrection();
	void ApplyCorrection(double g = 1.0);

	// Set the BPM data filter for this segment.
	// A NULL pointer indicates no filter is to
	// be applied.
	void SetBPMDataFilter(BPMDataFilter* f);

private:

	DFS_Segment itsSegment;
	BPMDataFilter* itsBPMdataFilter;

	ROChannelArray bpms;
	RWChannelArray correctors;
	RealVector cData;
	RealVector cCorr;
	TLAS::SVDMatrix<double>* svd;
	std::vector<RealVector> refdata;
};

#endif
