/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 3 - LatticeManipulation																	//
//																										//
//	The following provides a simple example of how to manually manipulate the constructed lattice		//
//	in order to simulate alignment errors/changes														//
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include units
#include "PhysicalUnits.h"

// Include MAD interface classes
#include "MADInterface.h"

// Include lattice manipulation classes
#include "StandardMultipoles.h"
#include "MagnetMover.h"

// Include closed orbit calculator
#include "LatticeFunctions.h"
#include "ClosedOrbit.h"

// Namespaces for convenience
using namespace PhysicalUnits;
using namespace ParticleTracking;

// Define accelerator parameters for convenience
#define beamenergy 5.0*GeV

int main()
{
	cout << "Locating MAD lattice information..." << endl;
	// Loop over possible build directories to locate MAD .tfs files
	string paths[] = {"../input/StorageRing.tfs", "input/StorageRing.tfs", "Tutorials/input/StorageRing.tfs"};
	string lattice_path;
	for (size_t i=0; i<3; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if (test_file)
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

	// Calculate beta and dispersion functions
	LatticeFunctionTable latticeFunctions(theModel,beamenergy);
	latticeFunctions.Calculate();

	// Write lattice functions to output file
	ofstream latticeFunctionLog("build/tutorial3a.out");
	latticeFunctions.PrintTable(latticeFunctionLog);

	// Manipulate lattice to simulation alignment errors etc
	vector<MagnetMover*> magnetMovers;
	theModel->ExtractTypedElements(magnetMovers);

	// Offset 20th magnet vertically by 100um
	magnetMovers[20]->SetY(100.0e-6);

	// Offset 40th magnet horizontally by 50um
	magnetMovers[40]->SetX(50.0e-6);

	// Extract a list of the quadrupoles
	vector<Quadrupole*> quadVec;
	theModel->ExtractTypedElements(quadVec);

	// Put a 20% gradient error on the 5th and 20th quadrupoles
	MultipoleField& field = quadVec[5]->GetField();
	Complex b1a = field.GetComponent(1);
	field.SetComponent(1, b1a.real() * 1.20, b1a.imag() * 1.20);

	MultipoleField& field2 = quadVec[20]->GetField();
	Complex b1b = field2.GetComponent(1);
	field2.SetComponent(1, b1b.real() * 1.10, b1b.imag() * 1.10);

	// Calculate beta and dispersion functions
	LatticeFunctionTable newLatticeFunctions(theModel,beamenergy);
	newLatticeFunctions.Calculate();

	// Write lattice functions to output file
	ofstream newLatticeFunctionLog("build/tutorial3b.out");
	newLatticeFunctions.PrintTable(newLatticeFunctionLog);

	delete theModel;

	cout << "Successful! Tutorial 3 Complete." << endl;
	cout << "Please run the corresponding python script to see a plot of manipulated lattice functions." << endl;

	return 0;
}
