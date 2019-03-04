/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "../tests.h"
#include <iostream>

#include "AcceleratorModelConstructor.h"
#include "Components.h"
#include "ParticleBunchTypes.h"
#include "ParticleTracker.h"
#include "CollimateParticleProcess.h"
#include "ApertureConfiguration.h"
#include "DetailedCollimationOutput.h"

/* Create a bunch of particles and track through a simple lattice with an aperture
 * output locations of losses to be checked by aperture_config_test.py
 */

using namespace std;

int main()
{
	AcceleratorModelConstructor ctor;

	AcceleratorComponent *drift;

	double s = 0;
	for(int i = 0; i < 6; i++)
	{
		double len = 1;
		drift = new Quadrupole("d" + to_string(i), len, 0);
		drift->SetComponentLatticePosition(s);
		ctor.AppendComponent(*drift);
		s += len;
	}

	AcceleratorModel* theModel = ctor.GetModel();

	double beam_energy = 7000.0;
	ProtonBunch* myBunch = new ProtonBunch(beam_energy, 1);

	size_t npart_r = 101, npart_t = 32;
	size_t npart = npart_r * npart_t;
	cout << "npart " << npart << endl;

	vector<Particle> pcoords;
	size_t id = 0;
	for(size_t i = 0; i < npart_r; i++)
	{
		double r = 1.0 / npart_r * i;
		for(size_t j = 0; j < npart_t; j++)
		{
			double t = 2 * M_PI / npart_t * j;
			Particle p(0);
			p.x() = r * sin(t);
			p.xp() = 0;
			p.y() = r * cos(t);
			p.yp() = 0;
			p.id() = id;
			id++;
			pcoords.push_back(p);
		}
	}

	AcceleratorModel::RingIterator ring = theModel->GetRing();
	ParticleTracker tracker(ring, myBunch);

	//#define COLL_AT_ENTRANCE 1
	//#define COLL_AT_CENTER 2
	//#define COLL_AT_EXIT 4
	CollimateParticleProcess* myCollimateProcess = new CollimateParticleProcess(2, 4);
	tracker.AddProcess(myCollimateProcess);

	DetailedCollimationOutput dco;
	dco.AddIdentifier("*");
	myCollimateProcess->SetCollimationOutput(&dco);

	ProtonBunch* test_bunch;
	vector<Particle> pcoords2;

	// Test with no aperture
	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy, 1, pcoords2);
	tracker.Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == npart);
	delete test_bunch;

	// Test with aperture from file
	ApertureConfiguration apc(find_data_file("aperture_config_test.tfs"));
	ofstream ap_log("act_aperture_conf.log");
	apc.SetLogFile(ap_log);
	apc.EnableLogging(1);
	apc.ConfigureElementApertures(theModel);
	cout << "aperture load finished" << endl << "start twiss" << endl;

	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy, 1, pcoords2);

	ofstream bunch_i("act_init_bunch.dat");
	test_bunch->Output(bunch_i, true);

	tracker.Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;

	ofstream losses("act_losses.dat");
	dco.Output(&losses);

	delete test_bunch;
	delete myBunch;
	delete theModel;
}
