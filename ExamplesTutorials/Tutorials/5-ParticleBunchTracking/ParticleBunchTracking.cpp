/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//																										//
//	Tutorial 5 - ParticleTracking																		//
//																										//
//	- Create a normally bistirbuted bunch of 10,000 particles and track around constructed lattice		//
//	- Record and plot BPM real coord centroid readings for 2 turns													//
//																										//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>

// Include units
#include "../ExamplesTutorials/Examples/Trajectory/BPMVectorBuffer.h"
#include "PhysicalUnits.h"

// Include MAD interface classes
#include "MADInterface.h"

// Include particle bunch classes
#include "BeamData.h"
#include "ParticleBunch.h"
#include "ParticleDistributionGenerator.h"
#include "ParticleTracker.h"

// Include BPM buffer

#define beamenergy 5.0 * GeV

using namespace PhysicalUnits;
using namespace ParticleTracking;

int main()
{
	size_t npart = (size_t) 1e4;

	int seed = 111;
	RandomNG::init(seed);

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

	AcceleratorModel* theModel = MADinput.ConstructModel();

	BeamData beam;
	beam.charge = 1e8;

	//TWISS beam parameters
	beam.beta_x = 0.5500005011 * meter;
	beam.beta_y = 0.5499999849 * meter;
	beam.alpha_x = -7.115569055e-7 * meter;
	beam.alpha_y = 1.797781918e-7 * meter;
	beam.emit_x = 5.026457122e-10 * meter;
	beam.emit_y = 5.026457122e-10 * meter;

	//Beam length.
	beam.sig_z = 75.5 * millimeter;

	//Relative energy spread of beam.
	beam.sig_dp = 0.000113;

	beam.p0 = beamenergy;

	ParticleDistributionGenerator* bunchDist = new NormalParticleDistributionGenerator();
	ParticleBunch* particleBunch = new ParticleBunch(npart, *bunchDist, beam);

	// Construct a ParticleTracker to perform the tracking
	ParticleTracker tracker(theModel->GetBeamline(), particleBunch);

	// Construct a BPMBuffer to record the bunch centroid at each BPM
	BPMVectorBuffer* bpmVecBuffer = new BPMVectorBuffer();
	BPM::SetDefaultBuffer(bpmVecBuffer);

	for(int turn = 0; turn < 2; turn++)
	{
		cout << "Tracking... turn: " << turn + 1 << endl;
		if(turn == 0)
		{
			tracker.Run();
		}
		else
		{
			tracker.Continue();
		}
		// Write the tracking results to a file
		ofstream bpmLog("build/tutorial5.out");
		vector<BPM::Data>& theBPMBuffer = bpmVecBuffer->BPMReading;
		for(vector<BPM::Data>::iterator bpm_iter = theBPMBuffer.begin(); bpm_iter != theBPMBuffer.end(); bpm_iter++)
		{
			bpmLog << std::setw(14) << (bpm_iter->x).value;
			bpmLog << std::setw(14) << (bpm_iter->x).error;
			bpmLog << std::setw(14) << (bpm_iter->y).value;
			bpmLog << std::setw(14) << (bpm_iter->y).error;
			bpmLog << endl;
		}
	}

	BPM::SetDefaultBuffer(nullptr);
	delete bpmVecBuffer;

	delete particleBunch;
	delete theModel;

	cout << "Successful! Tutorial 5 Complete." << endl;
	cout << "Please run the corresponding python script to see a plot of BPM real space coord readings." << endl;

	return 0;
}
