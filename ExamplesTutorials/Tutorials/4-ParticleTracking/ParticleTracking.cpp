/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 4 - ParticleTracking																		//
//																										//
//	- Create a bunch of 3 particles 1mm apart and track around constructed lattice						//
//	- Record and plot phase space coord variation after each turn for 100 turn							//
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>

// Include units
#include "PhysicalUnits.h"

// Include MAD interface classes
#include "MADInterface.h"

// Include lattice manipulation classes
#include "MagnetMover.h"

// Include closed orbit calculator
#include "ClosedOrbit.h"

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

	// Construct a bunch of 3 particles at 1mm horizontal intervals
	PSvector p(0);
	ParticleBunch* theParticles = new ParticleBunch(beamenergy);

	for(int xi = 1; xi <= 2; xi++)
	{
		p.x() = xi * 0.001;
		theParticles->AddParticle(p);
	}

	// Construct a ParticleTracker to perform the tracking
	ParticleTracker tracker(theModel->GetBeamline(), theParticles);

	// Track and record phase-space coords after each turn for 100 turns
	ofstream trackingLog("build/tutorial4.out");
	bool show_header = true;
	for(int turn = 0; turn < 100; turn++)
	{
		cout << "Tracking... turn: " << turn + 1 << endl;
		if(turn == 0)
		{
			tracker.Run();
		}
		else
		{
			show_header = false;
			tracker.Continue();
		}
		ParticleBunch& tracked = tracker.GetTrackedBunch();
		tracked.Output(trackingLog, show_header);
	}

	delete theParticles;
	delete theModel;

	cout << "Successful! Tutorial 4 Complete." << endl;
	cout << "Please run the corresponding python script to see a plot of tracked particles phase space coords." << endl;

	return 0;
}
