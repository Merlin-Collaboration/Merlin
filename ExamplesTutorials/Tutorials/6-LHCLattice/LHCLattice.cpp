/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 6 - The LHC																				//
//																										//
//	- Construct the entire LHC lattice, including aperture and collimator info in Merlin++				//
//	- Plot LHC lattice functions																		//
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include units and random number generator
#include "PhysicalUnits.h"
#include "RandomNG.h"

// Include MAD interface classes
#include "MADInterface.h"

// Include lattice function calculator
#include "LatticeFunctions.h"

// Include material, collimator and aperture config classes
#include "MaterialDatabase.h"
#include "CollimatorDatabase.h"
#include "ApertureConfiguration.h"

// Include closed orbit
#include "ClosedOrbit.h"

// Namespaces for convenience
using namespace PhysicalUnits;
using namespace ParticleTracking;

int main()
{
	//Define beam energy
	double beamenergy = 7000*GeV;

	// Import and construct LHC lattice
	cout << "Loading MAD lattice file..." << endl;
	MADInterface* MADinput = new MADInterface("ExamplesTutorials/Tutorials/input/LHC.tfs", beamenergy);
	MADinput->TreatTypeAsDrift("RFCAVITY");
	AcceleratorModel* theModel = MADinput->ConstructModel();

	// Import and define component aperture information
	cout << "Loading aperture information..." << endl;
	ApertureConfiguration* apertures = new ApertureConfiguration("ExamplesTutorials/Tutorials/input/LHCbeam1apertureinfo.tfs");
	apertures->ConfigureElementApertures(theModel);

	// Calculate closed orbit
	ClosedOrbit theClosedOrbit(theModel,beamenergy);
	Particle particle(0);
	theClosedOrbit.FindClosedOrbit(particle);

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

	// Define material database
	cout << "Loading materials database..." << endl;
	MaterialDatabase* material_db = new MaterialDatabase();

	// Import and define collimator information
	cout << "Loading collimators database..." << endl;
	CollimatorDatabase* collimator_db = new CollimatorDatabase("ExamplesTutorials/Tutorials/input/LHCcollimatorinfo.dat", material_db, true);

	// Write lattice functions to output file
	cout << "Writing lattice functions..." << endl;
	ofstream latticeFunctionLog("build/tutorial6.out");
	latticefunctions->PrintTable(latticeFunctionLog);

	delete collimator_db;
	delete material_db;
	delete apertures;
	delete theModel;
	delete MADinput;

	cout << "Successful! Tutorial 6 Complete." << endl;
	cout << "Please run the corresponding python script to see LHC lattice functions." << endl;

	return 0;
}
