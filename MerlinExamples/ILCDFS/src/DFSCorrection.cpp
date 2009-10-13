/////////////////////////////////////////////////////////////////////////
// Class DFSCorrection implementation
// 
// ILCDFS Application Code 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/19 10:19:06 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include "DFSCorrection.h"
#include "EnergyAdjustmentPolicy.h"
//#include "BPMDataFilter.h"
#include "TLAS/TLASimp.h"
#include "ResponseMatrixGenerator.h"
#include "ILCDFS_IO.h"

EnergyAdjustmentPolicy* DFSCorrection::theEnergyAdjustmentPolicy =0;
Accelerator* DFSCorrection::theReferenceModel =0;
Accelerator* DFSCorrection::theSimulationModel =0;

double DFSCorrection::w_abs = 1;
double DFSCorrection::w_diff =1;

DFSCorrection::DFSCorrection(const DFS_Segment &aSegment, Accelerator::Plane xy)
: itsSegment(aSegment), itsBPMdataFilter(0), bpms(), correctors(),svd(0),refdata()
{
	// First we use theReferenceModel to calculate the design 
	// response model matrix which will be used to calculate
	// future corrections on simulated data using theSimulationModel.
	dfs_trace(dfs_trace::level_1)<<"constructing DFS correction for segment "<<itsSegment<<' ';

	switch(xy) {
		case Accelerator::x_only: dfs_trace(dfs_trace::level_1)<<"X plane only"; break;
		case Accelerator::y_only: dfs_trace(dfs_trace::level_1)<<"Y plane only"; break;
		case Accelerator::x_and_y: dfs_trace(dfs_trace::level_1)<<"X and Y planes"; break;
	}
	dfs_trace(dfs_trace::level_1)<<endl;

	theReferenceModel->SetActiveBeamlineSegment(itsSegment);
	theEnergyAdjustmentPolicy->SetActiveBeamlineSegment(itsSegment);

	size_t nbpms = theReferenceModel->GetMonitorChannels(xy,bpms);
	size_t ncors = theReferenceModel->GetCorrectorChannels(xy,correctors);

	dfs_trace(dfs_trace::level_2)<<"No. of BPMs: "<<nbpms<<endl;
	dfs_trace(dfs_trace::level_2)<<"No. of correctors: "<<ncors<<endl;

	size_t nstates = theEnergyAdjustmentPolicy->GetNumEnergyStates();
	dfs_trace(dfs_trace::level_2)<<nstates<<"energy states"<<endl;
	
	RealMatrix M0(nbpms,ncors,0.0);
	RealMatrix Mi(nbpms,ncors,0.0);
	RealMatrix M(nstates*nbpms,ncors,0.0);
	refdata.reserve(nstates);

	dfs_trace(dfs_trace::level_1)<<"Constructing DFS model matrix"<<endl;

	ResponseMatrixGenerator rmg(theReferenceModel,bpms,correctors);

	for(size_t nes=0; nes<nstates; nes++) {

		dfs_trace(dfs_trace::level_2)<<"Constructing response matrix for state "<<nes<<endl;
		SubMatrix<double> m = M(Range(nes*nbpms,(nes+1)*nbpms-1),Range(0,ncors-1));		
		theEnergyAdjustmentPolicy->SetEnergyState(nes);

		rmg.Generate(nes);
		refdata.push_back(rmg.GetReference());

		if(nes==0){
			m = M0 = rmg.GetMatrix(); // constraint on absolute orbit
		}
		else {
			// For the off-energy states we take the 
			// difference matrix
			m = rmg.GetMatrix() - M0;
			refdata[nes] -= refdata[0]; // for off-energy states we need the difference orbit
		}
	}
	theEnergyAdjustmentPolicy->Restore();

	// Set up weight vector for SVD
	RealVector w(nstates*nbpms);
	w[Range(0,nbpms-1)] = w_abs;
	w[Range(nbpms,nstates*nbpms-1)] = w_diff;

	// Now construct the SVD of the complete model matrix
	dfs_trace(dfs_trace::level_2)<<"Performing SVD..."<<flush;
	svd = new TLAS::SVDMatrix<double>(M,w);
	dfs_trace(dfs_trace::level_2)<<"succesful"<<endl;

	// Finally we set set up the necessary channels 
	// for the actual simulated correction
	theSimulationModel->SetActiveBeamlineSegment(itsSegment);
	theSimulationModel->GetMonitorChannels(xy,bpms);
	theSimulationModel->GetCorrectorChannels(xy,correctors);
	cData.redim(nstates*nbpms);
	cCorr.redim(ncors);
};

void DFSCorrection::Initialise()
{
	dfs_trace(dfs_trace::level_2)<<"---- CORRECTING SEGMENT "<<itsSegment<<" --------- "<<endl;
	theSimulationModel->SetActiveBeamlineSegment(itsSegment);
	theEnergyAdjustmentPolicy->SetActiveBeamlineSegment(itsSegment);
}

void DFSCorrection::RecordTrajectories()
{
	dfs_trace(dfs_trace::level_2)<<"Recording data for segment "<<itsSegment<<endl;
	size_t nbpms = bpms.Size();

	// Simulated data acquisition
	size_t nstates = theEnergyAdjustmentPolicy->GetNumEnergyStates();
	RealVector xy0(bpms.Size());
	for(size_t nes=0; nes<nstates; nes++) {
		theEnergyAdjustmentPolicy->SetEnergyState(nes);
		theSimulationModel->TrackBeam(nes);
		if(nes==0) {
			xy0 = bpms;
			// difference to design absolute orbit
			cData(Range(0,nbpms-1)) = xy0-refdata[0]; 
		}
		else
			// difference from design 'difference orbit' for off-energy state
			cData(Range(nes*nbpms,(nes+1)*nbpms-1)) = (static_cast<RealVector>(bpms)-xy0)-refdata[nes];
	}
	theEnergyAdjustmentPolicy->Restore();

	if(dfs_trace::verbosity>=dfs_trace::level_3) {
		double rms = sqrt(cData*cData);
		dfs_trace(dfs_trace::level_3)<<"segment data rms = "<<rms<<endl;
	}
}

void DFSCorrection::CalculateCorrection()
{
	cCorr = (*svd)(cData);
}

void DFSCorrection::ApplyCorrection(double g)
{
	dfs_trace(dfs_trace::level_2)<<"Applying "<<100.0*g<<"% of correction"<<endl;
	RealVector x = -g*cCorr;
	correctors.IncrementAll(x);
	if(dfs_trace::verbosity>=dfs_trace::level_3) {
		RecordTrajectories();
	}

}

