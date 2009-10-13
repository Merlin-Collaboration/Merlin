/////////////////////////////////////////////////////////////////////////
// Class DFSApp implementation
// Root control class for applying DFS to a beamline
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

#include "DFSApp.h"
#include "DFSCorrection.h"
#include "EnergyAdjustmentPolicy.h"
#include "ILCDFS_IO.h"
#include "AcceleratorWithErrors.h"
#include "BeamDynamicsModel.h"

DFSApp::DFSApp()
{
	// should set-up defaults here
}

DFSApp::~DFSApp()
{
	// nothing to do
}

void DFSApp::SetReferenceModel(AcceleratorModel* accm, BeamData* beam0)
{
	Accelerator* acc = new Accelerator("REFERENCE MODEL",accm,beam0);
	DFSCorrection::theReferenceModel = acc;
}

void DFSApp::SetSimulationModel(AcceleratorModel* accm, BeamData* beam0)
{
	AcceleratorWithErrors* acc = new AcceleratorWithErrors("SIMULATION MODEL",accm,beam0);
	DFSCorrection::theSimulationModel = acc;
}

void DFSApp::SetReferenceModelBeamDynamics(BeamDynamicsModel* bdm)
{
	DFSCorrection::theReferenceModel->SetBeamDynamicsModel(bdm);
}

void DFSApp::SetSimulationModelBeamDynamicsForDFS(BeamDynamicsModel* bdm)
{
	dfsTracker = bdm;
}

void DFSApp::SetSimulationModelBeamDynamicsForEmittance(BeamDynamicsModel* bdm)
{
	emittanceTracker = bdm;
}

void DFSApp::SetEnergyAdjustmentPolicy(EnergyAdjustmentPolicy* epol)
{
	DFSCorrection::theEnergyAdjustmentPolicy = epol;
}

void DFSApp::SetSegmentation(size_t bBpmsPerSeg, size_t overlap)
{
	nbs = bBpmsPerSeg;
	nbo = overlap;
}

void DFSApp::SetIteration(size_t n)
{
	nit = n;
}

void DFSApp::SetCorrectionGain(double d)
{
	gain = d;
}

void DFSApp::SetCorrectionPlane(Accelerator::Plane p)
{
	pxy = p;
}

void DFSApp::Initialise()
{
	Accelerator* refacc = DFSCorrection::theReferenceModel;
	EnergyAdjustmentPolicy* engyPolicy = DFSCorrection::theEnergyAdjustmentPolicy;

	dfs_trace(dfs_trace::level_1)<<"Initialising DFS"<<endl;

	ReferenceParticleArray refpArray;
	KlystronArray klysArray;

	refacc->InitialiseTracking(engyPolicy->GetNumEnergyStates(),refpArray);
	refacc->GetKlystrons(klysArray);
	
	engyPolicy->SetKlystrons(klysArray);
	engyPolicy->SetReferenceParticles(refpArray);
	engyPolicy->Initialise();

	refacc->AllowIncrementalTracking(ref_useIT && engyPolicy->SupportsIncrementalTracking());

	DFS_Segment s = refacc->GetBeamlineRange();
	dfs_trace(dfs_trace::level_2)<<"beamline range: "<<s<<endl;

	IntegerArray ibpm;
	refacc->GetBeamlineIndecies("BPM.*",ibpm);

	dfs_trace(dfs_trace::level_1)<<ibpm.size()<<" BPMs in accelerator"<<endl;
	
	for(size_t n=0; n<ibpm.size(); n+=nbs-nbo) {
		DFS_Segment sn;
		sn.first = n==0 ? s.first : ibpm[n]+1;
		sn.second = n+nbs<ibpm.size() ? ibpm[n+nbs] : s.second;
		cachedDFS.push_back(new DFSCorrection(sn,pxy));
		if(sn.second == s.second)
			break;
	}
	dfs_trace(dfs_trace::level_1)<<cachedDFS.size()<<" segments initialised"<<endl;
}

void DFSApp::Apply()
{
	Accelerator* simacc = DFSCorrection::theSimulationModel;
	EnergyAdjustmentPolicy* engyPolicy = DFSCorrection::theEnergyAdjustmentPolicy;

	simacc->SetBeamDynamicsModel(dfsTracker);

	ReferenceParticleArray refpArray;
	KlystronArray klysArray;

	simacc->InitialiseTracking(engyPolicy->GetNumEnergyStates(),refpArray);
	simacc->GetKlystrons(klysArray);
	
	engyPolicy->SetKlystrons(klysArray);
	engyPolicy->SetReferenceParticles(refpArray);
	engyPolicy->Initialise();

	simacc->AllowIncrementalTracking(sim_useIT && engyPolicy->SupportsIncrementalTracking());

	for(list<DFSCorrection*>::iterator dfsi = cachedDFS.begin(); dfsi!=cachedDFS.end(); dfsi++) {
		(*dfsi)->Initialise();
		for(size_t n = 0; n<nit; n++) {
			dfs_trace(dfs_trace::level_2)<<"iteration "<<(n+1)<<'/'<<nit<<endl;
			(*dfsi)->RecordTrajectories();
			(*dfsi)->CalculateCorrection();
			(*dfsi)->ApplyCorrection(gain);
		}
	}

	engyPolicy->Restore();
	simacc->SetBeamDynamicsModel(emittanceTracker);
	simacc->TrackNewBunchThroughModel();
}

void DFSApp::ClearCachedDFS()
{
	for(list<DFSCorrection*>::iterator i = cachedDFS.begin(); i!=cachedDFS.end(); i++)
		delete *i;
	cachedDFS.clear();
}

AcceleratorWithErrors* DFSApp::GetSimulationModel()
{
	return static_cast<AcceleratorWithErrors*>(DFSCorrection::theSimulationModel);
}

Accelerator* DFSApp::GetReferenceModel()
{
	return DFSCorrection::theReferenceModel;
}

void DFSApp::SetOutput(SimulationOutput* dfso)
{
	emittanceTracker->SetOutput(dfso);
}

void DFSApp::SetFitWeights(double w_abs, double w_diff)
{
	DFSCorrection::w_abs = w_abs;
	DFSCorrection::w_diff = w_diff;
}

void DFSApp::UseIncrementalTracking(bool for_ref, bool for_sim)
{
	ref_useIT = for_ref;
	sim_useIT = for_sim;
}
