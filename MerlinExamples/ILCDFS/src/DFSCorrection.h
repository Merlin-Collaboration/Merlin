/////////////////////////////////////////////////////////////////////////
// Class DFSCorrection
// Applies the DFS correction to a segment of accelerator.
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/19 10:19:06 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_DFSCorrection
#define _h_DFSCorrection 1

#include "Channels/Channels.h"
#include "TLAS/LinearAlgebra.h"
#include "CommonDataStructures.h"
#include "Accelerator.h"

class EnergyAdjustmentPolicy;
class BPMDataFilter;
class SVDMatrix;

class DFSCorrection {
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
	void ApplyCorrection(double g=1.0);
	
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

