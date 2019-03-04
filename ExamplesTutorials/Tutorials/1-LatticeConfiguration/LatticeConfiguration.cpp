/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 1 - LatticeConfig																			//
//																										//
//	- The following provides a simple example of how to use the Merlin++ in-house lattice constructor	//
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Include i/o file stream
#include <fstream>

// Include units and constants
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

// Include model constructor headers
#include "AcceleratorModelConstructor.h"
#include "StandardMultipoles.h"
#include "SectorBend.h"

// Include closed orbit and lattice function calculators
#include "ClosedOrbit.h"
#include "LatticeFunctions.h"

// Namespaces for convenience
using namespace std;
using namespace PhysicalUnits;
using namespace PhysicalConstants;

// Define accelerator parameters for convenience
#define circum 1000.0
#define ncell 20.0
#define lcell 50.0
#define lquad 3
#define ldipole 5
#define dipolepercell 4
#define beamenergy 20*GeV
#define rigid beamenergy/eV/SpeedOfLight
#define curv 2*pi/(dipolepercell*ncell)

int main() {

	// Instantiate AcceleratorModelConstructor
	AcceleratorModelConstructor latticeConstructor;
	latticeConstructor.NewModel();

	cout << "Constructing Lattice..." << endl;

	// Period FODO lattice storage ring - loops over lattice cell ncell times
	//
	// WARNING: This is a basic and unstable FODO lattice without an RF component.
	// To resolve this issue for simple lattice analysis, we force longitudinal stability
	//
	for (int n=1;n<(ncell+1);++n) {
		latticeConstructor.AppendComponent(new Quadrupole("QF",lquad,0.0098*rigid), n==1 ? 0 : 0.15*lcell-ldipole);
		latticeConstructor.AppendComponent(new SectorBend("MB",ldipole,curv,rigid*curv), 0.15*lcell-lquad);
		latticeConstructor.AppendComponent(new SectorBend("MB",ldipole,curv,rigid*curv), 0.2*lcell-ldipole);
		latticeConstructor.AppendComponent(new Quadrupole("QD",lquad,-0.0098*rigid), 0.15*lcell-ldipole);
		latticeConstructor.AppendComponent(new SectorBend("MB",ldipole,curv,rigid*curv), 0.15*lcell-lquad);
		latticeConstructor.AppendComponent(new SectorBend("MB",ldipole,curv,rigid*curv), 0.2*lcell-ldipole);
	}
	latticeConstructor.AppendDrift(0.15*lcell-ldipole);
	AcceleratorModel* lattice = latticeConstructor.GetModel();

	cout << "Finding Closed Orbit..." << endl;

	// Find the closed orbit in the ring
	ClosedOrbit theClosedOrbit(lattice,beamenergy);
	Particle particle(0);
	theClosedOrbit.FindClosedOrbit(particle);

	// Calculate beta and dispersion functions
	LatticeFunctionTable* latticeFunctions = new LatticeFunctionTable(lattice,beamenergy);
	latticeFunctions->SetForceLongitudinalStability(true);
	latticeFunctions->Calculate();

	// Write lattice functions to output file
	ofstream latticeFunctionLog("build/tutorial1.out");
	latticeFunctions->PrintTable(latticeFunctionLog);

	delete lattice;

	cout << "Successful! Tutorial 1 Complete." << endl;

	return 0;
}


