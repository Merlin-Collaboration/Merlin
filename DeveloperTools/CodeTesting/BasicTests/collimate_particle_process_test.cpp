/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>

#include "AcceleratorModelConstructor.h"
#include "Components.h"
#include "PhysicalUnits.h"
#include "ParticleBunchTypes.h"
#include "ParticleTracker.h"
#include "Aperture.h"
#include "CollimateParticleProcess.h"

/* Create a bunch of particle, and check that the correct number survive various sized apertures
 *
 */

using namespace std;
using namespace PhysicalUnits;

int main(int argc, char* argv[])
{
	AcceleratorModelConstructor* ctor = new AcceleratorModelConstructor();
	ctor->NewModel();

	AcceleratorComponent *drift = new Drift("d1", 1 * meter);

	ApertureFactory factory;
	DataTable dt;
	dt.AddColumn("S", 'd');
	dt.AddColumn("APER_1", 'd');
	dt.AddColumn("APER_2", 'd');
	dt.AddColumn("APER_3", 'd');
	dt.AddColumn("APER_4", 'd');
	dt.AddColumn("APERTYPE", 's');
	dt.AddRow(0.0, 21 * millimeter, 100 * millimeter, 0.0, 0.0, "RECTELLIPSE");
	DataTableRowIterator itr = dt.begin();
	Aperture* rect_app = factory.GetInstance(*itr);
	ctor->AppendComponent(*drift);

	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;

	double beam_energy = 7000.0;
	ProtonBunch* myBunch = new ProtonBunch(beam_energy, 1);

	double coords[][4] =
	{
		{0, 0, 0, 0},
		{10, 0, 10, 0},
		{-20, 0, 20, 0},
		{10, 0, 0, 0},
		{20, 0, 0, 0},
		{0, 0, 10, 0},
		{0, 0, 20, 0},
		{0, 0, -30, 0},
	};
	size_t npart = sizeof(coords) / sizeof(coords[0]);
	cout << "npart " << npart << endl;

	vector<Particle> pcoords;
	for(size_t i = 0; i < npart; i++)
	{
		Particle p(0);
		p.x() = coords[i][0] * millimeter;
		p.xp() = coords[i][1];
		p.y() = coords[i][2] * millimeter;
		p.yp() = coords[i][3];
		pcoords.push_back(p);
		//myBunch->AddParticle(p);
	}

	AcceleratorModel::RingIterator ring = theModel->GetRing();
	ParticleTracker* tracker = new ParticleTracker(ring, myBunch);

	//#define COLL_AT_ENTRANCE 1
	//#define COLL_AT_CENTER 2
	//#define COLL_AT_EXIT 4
	CollimateParticleProcess* myCollimateProcess = new CollimateParticleProcess(2, 4);
	tracker->AddProcess(myCollimateProcess);

	drift->SetAperture(rect_app);

	ProtonBunch* test_bunch;
	vector<Particle> pcoords2;
	//
	// Test a few rectangular appertures with tracking
	//
	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy, 1, pcoords2);
	dt.AddRow(0.0, 35 * millimeter, 35 * millimeter, 0.0, 0.0, "RECTELLIPSE");
	itr++;
	Aperture* rect_app1 = factory.GetInstance(*itr);
	drift->SetAperture(rect_app1);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 8);
	delete test_bunch;

	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy, 1, pcoords2);
	dt.AddRow(0.0, 35 * millimeter, 25 * millimeter, 0.0, 0.0, "RECTELLIPSE");
	itr++;
	Aperture* rect_app2 = factory.GetInstance(*itr);
	drift->SetAperture(rect_app2);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 7);
	delete test_bunch;

	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy, 1, pcoords2);
	dt.AddRow(0.0, 25 * millimeter, 35 * millimeter, 0.0, 0.0, "RECTELLIPSE");
	itr++;
	Aperture* rect_app3 = factory.GetInstance(*itr);
	drift->SetAperture(rect_app3);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 8);
	delete test_bunch;

	pcoords2 = pcoords;
	test_bunch = new ProtonBunch(beam_energy, 1, pcoords2);
	dt.AddRow(0.0, 15 * millimeter, 15 * millimeter, 0.0, 0.0, "RECTELLIPSE");
	itr++;
	Aperture* rect_app4 = factory.GetInstance(*itr);
	drift->SetAperture(rect_app4);
	tracker->Track(test_bunch);
	cout << "Particle number: " << test_bunch->size() << endl;
	assert(test_bunch->size() == 4);
	delete test_bunch;

	delete rect_app;
	delete rect_app2;
	delete rect_app3;
	delete rect_app4;

	delete myBunch;
	delete tracker;
	delete theModel;

	cout << "test successful" << endl;

}
