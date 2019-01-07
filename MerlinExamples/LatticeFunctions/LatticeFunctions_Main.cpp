/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "PhysicalUnits.h"
#include "MADInterface.h"
#include "StandardMultipoles.h"
#include "MagnetMover.h"
#include "ClosedOrbit.h"
#include "LatticeFunctions.h"

#define BEAMENERGY 20 * GeV

typedef vector<MagnetMover*> MagnetMoverList;
typedef vector<Quadrupole*> QuadList;

using namespace PhysicalUnits;

int main()
{
	// Construct the AcceleratorModel
	// from a lattice file produced by MAD
	string paths[] = {"../lattices/twiss.out", "lattices/twiss.out", "MerlinExamples/lattices/twiss.out"};

	//twiss.7.0tev.b1_new.tfs
	string lattice_path;
	for(size_t i = 0; i < 3; i++)
	{
		ifstream test_file;
		test_file.open(paths[i].c_str());
		if(test_file)
		{
			lattice_path = paths[i];
			break;
		}
	}
	cout << "Lattice " << lattice_path << endl;

	MADInterface madi(lattice_path, BEAMENERGY);

	ofstream madlog("mad.log");
	madi.SetLogFile(madlog);
	madi.SetLoggingOn();

	AcceleratorModel* theModel = madi.ConstructModel();

	ClosedOrbit theClosedOrbit(theModel, BEAMENERGY);
	Particle co(0);
	theClosedOrbit.FindClosedOrbit(co);

//
//	// Extract a list of magnet movers from the AcceleratorModel
//	// and translate the 20th mover 20 microns vertically.
////	MagnetMoverList magnetMovers;
////	theModel->ExtractTypedElements(magnetMovers);
////	magnetMovers[20]->SetY(20.0e-6);
//
//
//	// Extract a list of the quadrupoles,
//	// and put a 5% gradient error on the 20th quadrupole.
////	QuadList quads;
////	theModel->ExtractTypedElements(quads);
////	MultipoleField& field = quads[20]->GetField();
////	Complex b1 = field.GetComponent(1);
////	field.SetComponent(1, b1.real()*1.05, b1.imag()*1.05);
//
//
//	// Find the lattice functions.
//	// Note the default is to find the closed orbit first,
//	// and use the transfer matrix around the closed orbit.
//	// Default functions are the closed orbit co-ordinates
//	// and the coupled-lattice equivalents of the
//	// Twiss alpha and beta functions.
	LatticeFunctionTable latticeFunctions(theModel, BEAMENERGY);

	latticeFunctions.UseOrbitFunctions();

//	// We add functions to give us the dispersion.
//	// Dx  = beta(1,6,3)/beta(6,6,3)
//	// Dpx = beta(2,6,3)/beta(6,6,3)
//	// etc.
//	latticeFunctions.AddFunction(1,6,3);
//	latticeFunctions.AddFunction(2,6,3);
//	latticeFunctions.AddFunction(3,6,3);
//	latticeFunctions.AddFunction(4,6,3);
//	latticeFunctions.AddFunction(6,6,3);
//
//	// Calculate the lattice functions.
	latticeFunctions.Calculate();
//
//	// Generate a file with the results.
//	// First column is the s position in the beamline.
//	ofstream latticeFunctionLog("LatticeFunctions.dat");
//	latticeFunctions.PrintTable(latticeFunctionLog);

	delete theModel;

	cout << "Finished!" << endl;
	return 0;
}
