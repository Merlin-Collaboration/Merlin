/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "BeamData.h"
#include "ParticleBunchConstructor.h"
#include "ParticleTracker.h"
#include "CollimatorWakeProcess.h"
#include "MADInterface.h"

#include "RandomNG.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

#include "SynchRadParticleProcess.h"
#include "CollimateParticleProcess.h"

#include "CollimatorDatabase.h"
#include "MaterialDatabase.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;

int main(int argc, char* argv[])
{

	MaterialDatabase* mat = new MaterialDatabase();

	CollimatorDatabase* collimator_db = new CollimatorDatabase("collimator.db", mat, 0);

	int seed = 0;
	if(argc >= 2)
	{
		seed = atoi(argv[1]);
	}
	cout << "Random Seed: " << seed << endl;

	//Initialise Random number generator
	RandomNG::init(seed);

	double offset;
	int modes;
	modes = 1;
	offset = 0;

	int npart = 1000;
	int nlaps = 10;

	/*********************************************************************
	 *
	 *
	 *	BEAM SETTINGS
	 *
	 *
	 *********************************************************************/

	//Create a beam
	BeamData mybeam;

	//Default values are 0.0

	//The charge of the particles in the beam.
	//  <0 for electrons, >0 for positrons/protons.
	mybeam.charge = 1.11e11;

	//TWISS beam parameters
	mybeam.beta_x = 0.5500005011 * meter;
	mybeam.beta_y = 0.5499999849 * meter;
	mybeam.alpha_x = -7.115569055e-7 * meter;
	mybeam.alpha_y = +1.797781918e-7 * meter;
	mybeam.emit_x = 5.026457122e-10 * meter;
	mybeam.emit_y = 5.026457122e-10 * meter;

	//Beam length.
	mybeam.sig_z = 75.5 * millimeter;

	//Relative energy spread of beam.
	mybeam.sig_dp = 0.000113;

	//Beam centroid
	mybeam.x0 = 0;
	mybeam.xp0 = 0;
	mybeam.y0 = 0;
	mybeam.yp0 = 0;
	mybeam.ct0 = 0.0;

	//Beam energy (momentum).
	mybeam.p0 = 7000 * GeV;
	//cout << "Energy = " << mybeam.p0 <<endl;

	//X-Y coupling
	mybeam.c_xy = 0.0;
	mybeam.c_xyp = 0.0;
	mybeam.c_xpy = 0.0;
	mybeam.c_xpyp = 0.0;

	//Dispersion
	mybeam.Dx = 0.0;
	mybeam.Dxp = 0.0;
	mybeam.Dy = 0.0;
	mybeam.Dyp = 0.0;

	//Check if beam parameters are ok.
	if(!mybeam.ok())
	{
		cerr << "Bad beam parameters: Check emittance and beta." << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Beam parameters ok." << endl;

	/*********************************************************************
	 *
	 *
	 *	BUNCH SETTINGS
	 *
	 *
	 *********************************************************************/

	ParticleBunch* myBunch = ParticleBunchConstructor(mybeam, npart, normalDistribution).ConstructParticleBunch();
	//ParticleBunch* myBunch = ParticleBunchConstructor(mybeam,npart,pencilDistribution).ConstructParticleBunch();

	//Output the initial input bunch
	ofstream *bunchbeforefile = new ofstream("Output/bunch0.dat");

	myBunch->Output(*bunchbeforefile);
	delete bunchbeforefile;

	/*********************************************************************
	 *
	 *
	 *	ACCELERATOR MODEL LOADING
	 *
	 *
	 *********************************************************************/

	//Load accelerator optics file.
	MADInterface* myMADinterface = new MADInterface("MerlinExamples/LHC/Input/lhc.tfs", 7000.0 * GeV);
	//MADInterface* myMADinterface = new MADInterface("Input/LHCB19.tfs",7000.0*GeV);
	cout << "Got MADInterface" << endl;

	//Enable/Disable apertures
	//myMADinterface->incApertures=true;

	//Set MADInterface log file
	ofstream MADLog("Output/MADlog.txt");
	myMADinterface->SetLogFile(MADLog);

	//Enable Logging
	myMADinterface->SetLoggingOn();

	//Set the Collimator db
	//myMADinterface->Set_Collimator_Database(collimator_db);

	//Build accelerator model
	AcceleratorModel* model = myMADinterface->ConstructModel();
	cout << "Built MADInterface" << endl;

	//Output the accelerator model component statistics
	ofstream myoutfile("Output/model.txt");
	model->ReportModelStatistics(myoutfile);

	/*********************************************************************
	 *
	 *
	 *	PARTICLE TRACKER
	 *
	 *
	 *********************************************************************/

	ParticleTracker* tracker = new ParticleTracker(model->GetBeamline(), myBunch);
	//ParticleTracker* tracker = new ParticleTracker(model->GetRing(),myBunch);

	/*********************************************************************
	 *
	 *
	 *	COLLIMATION SETTINGS
	 *
	 *
	 *********************************************************************/

	//Output stream for collimator losses
	ofstream* myout = new ofstream("Output/loss.txt");
	if(!myout->good())
	{
		std::cerr << "Could not open collimation loss file" << std::endl;
		exit(EXIT_FAILURE);
	}

	//New Collimation process
	CollimateParticleProcess* myCollimateProcess = new CollimateParticleProcess(2, 7, myout);

	//Enable scattering
	myCollimateProcess->ScatterAtCollimator(true);

	//Create individual loss files
	myCollimateProcess->CreateParticleLossFiles(true, "lostlist");

	// Sets maximum allowed loss percentage at a single collimator.
	myCollimateProcess->SetLossThreshold(101.0);

	//Output stream for collimator process log
	ofstream* myCollimationLog = new ofstream("Output/collimator_log.txt");

	//sets process log stream, NULL to disable.
	myCollimateProcess->SetLogStream(NULL);

	//Add Collimation process to the tracker.
	tracker->AddProcess(myCollimateProcess);

	/*********************************************************************
	 *
	 *
	 *	WAKEFIELD SETTINGS
	 *
	 *
	 *********************************************************************/

	// apply the resistive wakefields
	// modes, priority, nbins, nsigma
	CollimatorWakeProcess* myWakeProcess = new CollimatorWakeProcess(modes, 1, 10, 3);

	//Enable the Wakefield process
	tracker->AddProcess(myWakeProcess);

	/*********************************************************************
	 *
	 *
	 *	SYNCHROTRON RADIATION SETTINGS
	 *
	 *
	 *********************************************************************/

	//SynchRadParticleProcess* mySynchRadParticleProcess = new SynchRadParticleProcess(10, 1);

	//Include radiation effects in Quadrupoles and Skew Quadrupoles. Bool switch.
	//mySynchRadParticleProcess->IncludeQuadRadiation (1);

	//Set photon generation type HBSpectrumGen or AWSpectrumGen
	//mySynchRadParticleProcess->SetPhotonGenerator(HBSpectrumGen);

	//Set number of steps though each component, default = 1
	//mySynchRadParticleProcess->SetNumComponentSteps(1);

	//Enable Process - currently will segfault
	//tracker ->AddProcess(mySynchRadParticleProcess);

	// Do the loop for nlaps times
	for(int iii = 1; iii <= nlaps; iii++)
	{
		ParticleArray testarray = myBunch->GetParticles();
		cout << "Lap " << iii << "\tParticle number: " << testarray.size() << endl;
		ostringstream lapprefix;
		lapprefix << "Output/loss_lap_" << setw(4) << setfill('0') << iii << "_";
		myCollimateProcess->CreateParticleLossFiles(true, lapprefix.str().c_str());
		tracker->Track(myBunch);
	}

	cout << "Laps done" << endl;
	ofstream *bunchafterfile = new ofstream("Output/bunch1.dat");
	myBunch->Output(*bunchafterfile);
	delete bunchafterfile;

	return EXIT_SUCCESS;
}

//mybeam.beta_x = 0.5495121695*meter;
//mybeam.beta_y = 0.5498820579*meter;
//mybeam.emit_x = 33.640*5.026457122e-10*meter;
//mybeam.emit_y = 33.64*5.026457122e-10*meter;
//mybeam.alpha_x = -0.0001721885021*meter;
//mybeam.alpha_y = -0.0004654580947*meter;

/*	double benergy;
    ofstream energy_out("Output/Energy.txt");
    energy_out.precision(100);
    cout.precision(100);
 */
//	simp << testarray[1].x() << "\t" << testarray[1].xp() << endl;
//		benergy = myBunch->AdjustRefMomentumToMean();
//		energy_out << iii << "\t" << benergy << endl;

//Cleaning up
//delete myWakeProcess;
//delete myCollimateProcess;
//delete mySynchRadParticleProcess;

//delete myMADinterface;

//delete myCollimationLog;

//Does not like this being deleted
//delete tracker;

//delete zero_particle;
//delete myout;
