/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

// Beam Delivery System Collimator Studies
// --------------------------------------------------------------------
//
// The application tracks a halo of particles with the specified extents,
// and performs hard-aperture cuts where apertures are defined. Particles
// outside of a given aperture in the model are deemed lost, and are output
// to a particle loss file of the form
//
//   data/component_type.component_label.n.loss
//
// where n is the occurance number of the component in the beamline.
// examples:
//       data/Drift.SPOIX.1.loss
//       data/Drift.SPOIY.3.loss
//
// class HaloTracker also gives an example of how to override the default
// particle tracking mechanism (see HaloTracker.cpp and QuadIntegrator.[h,cpp]

#include "RandomNG.h"

// For this example, we incapsulate the halo construction
// and tracking in a utility class HaloTracker.
#include "HaloTracker.h"
#include <iostream>

#include "XTFFInterface.h"

using namespace std;

// --------------------------------------------------------------------
// global functions forward declarations
// --------------------------------------------------------------------

pair<AcceleratorModel*, BeamData*> ConstructModel(const string& fname);

// --------------------------------------------------------------------
// Main function
// --------------------------------------------------------------------

int main()
{
	// Initialise the random number generated FIRST!
	RandomNG::init();

	// Construct the BDS beamline model
	string paths[] = {"../lattices/tesla_bds_v8.05.optics", "lattices/tesla_bds_v8.05.optics",
					  "MerlinExamples/lattices/tesla_bds_v8.05.optics"};

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
	pair<AcceleratorModel*, BeamData*> mb = ConstructModel(lattice_path);

	AcceleratorModel* model = mb.first;
	BeamData* beam = mb.second;

	// Construct the halo tracking object.
	HaloTracker ht(model->GetBeamline(), *beam);

	// set halo limits to ±5 sigma and ±5% in dp/p
	ht.SetHaloLimitsN(5, 5, 5, 5, 0.05);

	// We want to hard collimate the halo at apertures,
	ht.collimate_halo = true;

	// Once tracking is complete, dump out
	// the surviving particles (data/particles.out.dat)
	ht.dump_particles = true;

	// Turn scattering on/off at Collimator elements
	ht.scatter_at_collimator = false;

	// Track halo with 1e4 particles
	ht.Run(10000);

	// clean up
	delete model;
	delete beam;

	return 0;
}
