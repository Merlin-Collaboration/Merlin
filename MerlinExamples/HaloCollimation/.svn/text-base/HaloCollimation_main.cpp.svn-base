// Beam Delivery System Spoiler Studies
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
//

#include "Random/RandomNG.h"

// For this example, we incapsulate the halo construction
// and tracking in a utility class HaloTracker.
#include "HaloTracker.h"

using namespace std;

// --------------------------------------------------------------------
// global functions forward declarations
// --------------------------------------------------------------------

pair<AcceleratorModel*,BeamData*> ConstructModel(const string& fname);

// --------------------------------------------------------------------
// Main function
// --------------------------------------------------------------------

int main()
{
    // Initialise the random number generated FIRST!
    RandomNG::init();

    // Construct the BDS beamline model
    pair<AcceleratorModel*,BeamData*> mb = ConstructModel("../lattices/tesla_bds_v8.05.optics");

    AcceleratorModel* model = mb.first;
    BeamData* beam = mb.second;
               
    // Construct the halo tracking object.
    HaloTracker ht(model->GetBeamline(),*beam);

    // set halo limits to ±5 sigma and ±5% in dp/p
    ht.SetHaloLimitsN(5,5,5,5,0.05);

    // We want to hard collimate the halo at apertures,
    ht.collimate_halo = true;

    // Once tracking is complete, dump out
    // the surviving particles (data/particles.out.dat)
    ht.dump_particles = true;

	// Turn scattering on/off at SPOILER elements
	ht.scatter_at_spoiler = false;

    // Track halo with 1e4 particles
    ht.Run(10000);
    
    // clean up
    delete model;
    delete beam;

    return 0;
}

