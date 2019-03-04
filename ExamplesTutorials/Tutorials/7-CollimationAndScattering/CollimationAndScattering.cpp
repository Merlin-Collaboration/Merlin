/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 7 - The LHC																				//
//																										//
//	- Construct the entire LHC lattice, including aperture and collimator info in Merlin++				//
//	- Define beam																						//
//	- Initialize collimation and scattering																//
//	- Record losses and plot loss map																	//
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include units and constants etc
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "RandomNG.h"
#include "NANCheckProcess.h"

// Include MAD interface classes
#include "MADInterface.h"

// Include lattice function calculator
#include "LatticeFunctions.h"
#include "Dispersion.h"

// Include particle bunch classes
#include "BeamData.h"
#include "ParticleBunch.h"
#include "ParticleBunchTypes.h"
#include "ParticleDistributionGenerator.h"
#include "HaloParticleDistributionGenerator.h"
#include "ParticleTracker.h"
#include "BunchFilter.h"

// Include collimator classes
#include "MaterialDatabase.h"
#include "CollimatorDatabase.h"
#include "CollimateProtonProcess.h"
#include "LossMapCollimationOutput.h"

// Include aperture classes
#include "ApertureConfiguration.h"
#include "CollimatorAperture.h"

// Include scattering classes
#include "ScatteringProcess.h"
#include "ScatteringModelsMerlin.h"

// Namespaces for convenience
using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

int main()
{
	// Define number of particles and turns
	size_t npart = (size_t) 1e4;
	size_t nturns = (size_t) 20;
	int seed = 111;
	RandomNG::init(seed);

	// Define main beam params
	double beamenergy = 7000*GeV;
	double beamcharge = 1.1e11;
	double normalized_emittance = 3.5e-6;
	double gamma = beamenergy / ProtonMassMeV / MeV;
	double beta = sqrt(1.0 - (1.0 / pow(gamma, 2)));
	double emittance = normalized_emittance / (gamma * beta);

	// Import and construct LHC lattice
	cout << "Loading MAD lattice file..." << endl;
	MADInterface* MADinput = new MADInterface("ExamplesTutorials/Tutorials/input/LHC.tfs", beamenergy);
	MADinput->TreatTypeAsDrift("RFCAVITY");
	AcceleratorModel* theModel = MADinput->ConstructModel();

	// Import and define component aperture information
	cout << "Loading aperture information..." << endl;
	ApertureConfiguration* apertures = new ApertureConfiguration("ExamplesTutorials/Tutorials/input/LHCbeam1apertureinfo.tfs");
	apertures->ConfigureElementApertures(theModel);

	// Define material database
	cout << "Loading materials database..." << endl;
	MaterialDatabase* material_db = new MaterialDatabase();

	// Import and define collimator information
	cout << "Loading collimators database..." << endl;
	CollimatorDatabase* collimator_db = new CollimatorDatabase("ExamplesTutorials/Tutorials/input/LHCcollimatorinfo.dat", material_db, true);

	// Instantiate and calculate lattice functions
	cout << "Calculating lattice functions..." << endl;
	LatticeFunctionTable* latticefunctions = new LatticeFunctionTable(theModel, beamenergy);

	// Dynamic lattice function convergence loop
	double bscale1 = 1e-22;
	while(true)
	{
		latticefunctions->ScaleBendPathLength(bscale1);
		latticefunctions->Calculate();
		if(!std::isnan(latticefunctions->Value(1, 1, 1, 0)))
		{
			break;
		}
		bscale1 *= 2;
	}

	// Calculate Dispersion
	string start_element = "TCP.C6L7.B1";
	int start_element_number = 8764;
	Dispersion* dispersion = new Dispersion(theModel, beamenergy);
	dispersion->FindDispersion(start_element_number);

	// Initialize collimator database
	collimator_db->MatchBeamEnvelope(false);
	collimator_db->EnableJawAlignmentErrors(false);
	collimator_db->SetJawPositionError(0.0 * nanometer);
	collimator_db->SetJawAngleError(0.0 * microradian);
	collimator_db->SelectImpactFactor(start_element, 1.0e-6);
	double impact = collimator_db->ConfigureCollimators(theModel, emittance, emittance, latticefunctions);

	// Fully initialize beam data using the above
	BeamData beamData;

	beamData.charge = beamcharge / npart;
	beamData.p0 = beamenergy;
	beamData.beta_x = latticefunctions->Value(1, 1, 1, start_element_number) * meter;
	beamData.beta_y = latticefunctions->Value(3, 3, 2, start_element_number) * meter;
	beamData.alpha_x = -latticefunctions->Value(1, 2, 1, start_element_number);
	beamData.alpha_y = -latticefunctions->Value(3, 4, 2, start_element_number);

	//Dispersion
	beamData.Dx = dispersion->Dx;
	beamData.Dy = dispersion->Dy;
	//mybeam.Dy=0;
	beamData.Dxp = dispersion->Dxp;
	beamData.Dyp = dispersion->Dyp;

	beamData.emit_x = impact * impact * emittance * meter;
	impact = 1;
	beamData.emit_y = impact * impact * emittance * meter;
	beamData.sig_z = 0.0;

	//Beam centroid
	beamData.x0 = latticefunctions->Value(1, 0, 0, start_element_number);
	beamData.xp0 = latticefunctions->Value(2, 0, 0, start_element_number);
	beamData.y0 = latticefunctions->Value(3, 0, 0, start_element_number);
	beamData.yp0 = latticefunctions->Value(4, 0, 0, start_element_number);
	beamData.ct0 = latticefunctions->Value(5, 0, 0, start_element_number);

	beamData.sig_dp = 0.0;

	//X-Y coupling
	beamData.c_xy = 0.0;
	beamData.c_xyp = 0.0;
	beamData.c_xpy = 0.0;
	beamData.c_xpyp = 0.0;

	// Initialize collimator aperture info
	vector<Collimator*> TCP;
	int siz = theModel->ExtractTypedElements(TCP, start_element);
	Aperture *ap = (TCP[0])->GetAperture();
	CollimatorAperture* CollimatorJaw = dynamic_cast<CollimatorAperture*>(ap);
	double h_offset = latticefunctions->Value(1, 0, 0, start_element_number);
	double JawPosition = CollimatorJaw->GetFullEntranceWidth() / 2.0;

	// Define horizontal bunch filter
	HorizontalHaloParticleBunchFilter* hFilter = new HorizontalHaloParticleBunchFilter();
	hFilter->SetHorizontalLimit(JawPosition);
	hFilter->SetHorizontalOrbit(h_offset);

	// Construct corresponding bunch
	ProtonBunch* particleBunch = new ProtonBunch(npart, HorizonalHalo2ParticleDistributionGenerator(), beamData, hFilter);
	particleBunch->SetMacroParticleCharge(beamData.charge);

	// Enable scattering physics
	particleBunch->EnableScatteringPhysics(ProtonBunch::Merlin);

	// Construct a ParticleTracker to perform tracking
	AcceleratorModel::RingIterator ring = theModel->GetRing(start_element_number);
	ParticleTracker* tracker = new ParticleTracker(ring, particleBunch);
	tracker->SetLogStream(cout);

	auto nancheck = new NANCheckProcess;
	nancheck->SetDetailed(0);
	nancheck->SetHaltNAN(1);
	tracker->AddProcess(nancheck);

	// Define collimation and scattering settings
	string filename = "build/tutorial7.out";

	LossMapCollimationOutput* lossOutput = new LossMapCollimationOutput(tencm);
	ScatteringModel* scatterModel = new ScatteringModelMerlin;

	CollimateProtonProcess* collimateProcess = new CollimateProtonProcess(2, 4);
	collimateProcess->SetScatteringModel(scatterModel);

	collimateProcess->ScatterAtCollimator(true);
	collimateProcess->SetLossThreshold(200.0);
	collimateProcess->SetOutputBinSize(0.1);
	collimateProcess->SetCollimationOutput(lossOutput);
	tracker->AddProcess(collimateProcess);

	// Run tracker
	for(size_t turn = 1; turn <= nturns; ++turn)
	{
		cout << "Turn " << turn << "\tParticle number: " << particleBunch->size() << endl;
		tracker->Track(particleBunch);
		if(particleBunch->size() <= 1)
		{
			break;
		}
	}

	lossOutput->Finalise();
	ofstream* col_output = new ofstream(filename);
	lossOutput->Output(col_output);
	col_output->flush();
	col_output->close();
	cout << "Output: " << filename << endl;

	cout << "npart: " << npart << endl;
	cout << "left: " << particleBunch->size() << endl;
	cout << "absorbed: " << npart - particleBunch->size() << endl;

	delete MADinput;
	delete theModel;
	delete apertures;
	delete material_db;
	delete collimator_db;
	delete latticefunctions;
	delete hFilter;
	delete dispersion;
	delete particleBunch;
	delete tracker;
	delete scatterModel;
	delete collimateProcess;
	delete col_output;

	cout << "Successful! Tutorial 7 Complete." << endl;
	cout << "Please run the corresponding python script to see tracked particle info and scattering loss maps." << endl;

	return 0;
}
