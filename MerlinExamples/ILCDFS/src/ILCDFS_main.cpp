/////////////////////////////////////////////////////////////////////////
// ILCDFS_main
// Main programme for ILC DFS simulation appliction
// 
// Based on the MERLIN class library
//
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2006/06/12 14:30:09 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#include "Random/RandomNG.h"
#include "DFSApp.h"
#include "ParticleTrackingModel.h"
#include "SMPTrackingModel.h"
#include "AcceleratorWithErrors.h"
#include "ConstantGradientAdjustment.h"
#include "ModelConstruction.h"
#include "ILCDFS_IO.h"
#include "NumericalUtils/PhysicalUnits.h"
#include "DFSOutput.h"
#include "OneToOneCorrection.h"

using namespace PhysicalUnits;

int main()
{
	RandomNG::init(555);
	DFSApp theDFSapp;

	// General DFS parameters
	// --------------------------------------------------------------------
	theDFSapp.SetCorrectionGain(1.0);
	theDFSapp.SetCorrectionPlane(Accelerator::y_only);
	theDFSapp.SetIteration(1.0);
	theDFSapp.SetSegmentation(40,20);
	theDFSapp.SetFitWeights(1/(360.0*micrometer),1/(sqrt(2.0)*5*micrometer));
	theDFSapp.UseIncrementalTracking(true,true);

//	theDFSapp.SetFitWeights(1,0); // 1-1 steering
//	theDFSapp.SetFitWeights(0,1); // dispersion correction

	// trace detail
	dfs_trace::verbosity = dfs_trace::level_1;
//	ofstream traceOS("trace.dat");
//	dfs_trace::set_trace_stream(traceOS);

	// Model construction
	// --------------------------------------------------------------------
	// Construct the simulation and reference
	// accelerator models.
	const string modelFile = "../lattices/ilc_linac_15_250.xtff";
	const bool constructCurvedLinac = true;

	pair<AcceleratorModel*,BeamData*> mb;
	mb = ConstructModel(modelFile,constructCurvedLinac);  
	theDFSapp.SetReferenceModel(mb.first,mb.second);
	mb = ConstructModel(modelFile,constructCurvedLinac);
	theDFSapp.SetSimulationModel(mb.first,mb.second);

	// Tracking
	// --------------------------------------------------------------------
	// Reference model uses single-ray tracing to simulate trajectory
	// (no wakefields)
	theDFSapp.SetReferenceModelBeamDynamics(new ParticleTrackingModel(1));

	// Simulation model uses single-ray tracing to simulate trajectory
	// during DFS application (no wakefields);
	theDFSapp.SetSimulationModelBeamDynamicsForDFS(new ParticleTrackingModel(1));

	// Simulation model uses SMPTracking for emittance estimation after 
	// DFS application.
	// 31 slices with 11 macro-particles per slice.
	SMPTrackingModel* smpTracker = new SMPTrackingModel(31,11);
	smpTracker->IncludeTransverseWakefield(true);
	theDFSapp.SetSimulationModelBeamDynamicsForEmittance(smpTracker);

	// Set-up the required output. For this example we will produce
	// a separate file for each seed containing the data at every BPM
	DFSOutput dfso;
	dfso.output_initial = true;
	dfso.output_final = true;
	dfso.output_all = false;
	dfso.AddIdentifier("BPM.*");

	// Accelerator (random) errors
	// --------------------------------------------------------------------
	AcceleratorWithErrors* accWE = theDFSapp.GetSimulationModel();

	accWE->TransverseErrors("ComponentFrame.Quadrupole.*",0,300*micrometer, true);
	accWE->RotationErrors("ComponentFrame.Quadrupole.*",0,0,300*microradian, false);
	accWE->TransverseErrors("ComponentFrame.BPM.*",0,200*micrometer, true);
	accWE->TransverseErrors("ComponentFrame.TWRFStructure.*",0,300*micrometer, true);
	accWE->RotationErrors("ComponentFrame.TWRFStructure.*",300*microradian,0,0, false);
	accWE->TransverseErrors("SequenceFrame.CM",0,200*micrometer, true); // Cryomodule
	accWE->MagnetScaleError("YCor.*",0.01);
	accWE->KlystronErrors("*",0.01,0.0);
	accWE->BPMresolution(5.0*micrometer);
	accWE->BPMlinearError(0.05);
	accWE->InitialBeamJitterInSigma(1.0,1.0);

	// Energy modification for DFS
	// --------------------------------------------------------------------
	ConstantGradientAdjustment* engyPolicy = new ConstantGradientAdjustment();
//	engyPolicy->AddEnergyState(0,-0.2); // state 1: -20% reduction in initial beam energy
//	engyPolicy->AddEnergyState(-0.2,0); // state 2: -20% reduction in linac gradient
	engyPolicy->AddEnergyState(-0.2,-0.2); // D. Schulte's approach.
	theDFSapp.SetEnergyAdjustmentPolicy(engyPolicy);

	// Simulation loop
	size_t nseed = 100;
	const string filePrefix = "test_results/seed";

	// If we are following the Earth's curvature, we must first one-to-one steer
	// the reference model before using it to calculate the reference trajectories
	// and response matrices
	if(constructCurvedLinac)
		OneToOneCorrection(theDFSapp.GetReferenceModel(),Accelerator::y_only);

	theDFSapp.Initialise(); // one-time initialisation for all seeds.
	theDFSapp.SetOutput(&dfso);
	dfs_trace(dfs_trace::level_1)<<"\n\n\n\nBeginning simulations with "<<nseed<<" seeds"<<endl;

	for(size_t ns = 0; ns<nseed; ns++) {
		dfs_trace(dfs_trace::level_1)<<"\nRunning seed #"<<ns;
		ostringstream oss;
		oss<<filePrefix<<'.'<<ns<<".dat";
		dfs_trace(dfs_trace::level_1)<<" file: "<<oss.str()<<endl;
		dfso.NewFile(oss.str());
		accWE->ApplyStaticErrors();
		theDFSapp.Apply();
	}

	return 0;
}

