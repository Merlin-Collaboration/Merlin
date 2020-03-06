/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 2 - LatticeConfigMAD																		//
//																										//
//	- The following provides a simple example of how to import a MAD lattice .tfs file into Merlin++	//
//	- A corresponding python script is provided to plot the calculated beta and dispersion functions    //
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include i/o file stream
#include <fstream>

// Include units
#include "PhysicalUnits.h"

// Include MAD interface classes
#include "MADInterface.h"

// Include closed orbit and lattice function calculators
#include "ClosedOrbit.h"
#include "LatticeFunctions.h"

// Namespaces for convenience
using namespace std;
using namespace PhysicalUnits;
using namespace ParticleTracking;

// Define accelerator parameters for convenience
#define beamenergy 5.0 * GeV

int main()
{
	cout << "Locating MAD lattice information..." << endl;
	// Loop over possible build directories to locate MAD .tfs files
	string paths[] = {"../input/StorageRing.tfs", "input/StorageRing.tfs", "Tutorials/input/StorageRing.tfs",
					  "ExamplesTutorials/Tutorials/input/StorageRing.tfs"};
	string lattice_path;
	for(size_t i = 0; i < 4; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if(test_file)
		{
			lattice_path = paths[i];
			break;
		}
	}
	cout << "MAD lattice information found at: " << lattice_path << endl;

	cout << "Constructing Lattice..." << endl;

	// Instantiate MADInterface with lattice file input
	MADInterface MADinput(lattice_path, beamenergy);

	// Construct Model
	AcceleratorModel* theModel = MADinput.ConstructModel();

	cout << "Finding Closed Orbit..." << endl;

	// Find the closed orbit in the ring
	ClosedOrbit theClosedOrbit(theModel, beamenergy);
	Particle particle(0);
	theClosedOrbit.FindClosedOrbit(particle);

	// Calculate beta and dispersion functions
	LatticeFunctionTable* latticeFunctions = new LatticeFunctionTable(theModel, beamenergy);
	//latticeFunctions->SetForceLongitudinalStability(true);
	latticeFunctions->Calculate();

	// Write lattice functions to output file
	ofstream latticeFunctionLog("build/tutorial2.out");
	latticeFunctions->PrintTable(latticeFunctionLog);

	delete theModel;

	cout << "Successful! Tutorial 2 Complete." << endl;
	cout << "Please run the corresponding python script to see a plot of the beta and dispersion functions." << endl;

	return 0;
}
