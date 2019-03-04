/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>
#include <fstream>

#include "AcceleratorModelConstructor.h"
#include "Components.h"
#include "PhysicalUnits.h"
#include "ParticleBunchTypes.h"
#include "ParticleTracker.h"
#include "HollowElectronLens.h"
#include "HollowELensProcess.h"
#include "NANCheckProcess.h"

/*
 * Test that the deterministic modes of the the Hollow Election Lens
 * agree with a reference file.
 *
 */

using namespace std;
using namespace PhysicalUnits;

int main()
{
	const double beam_energy = 7000.0;
	const size_t npart = 400;

	// Test bunch with particles along x, y and diagonals
	vector<Particle> pcoords;
	for(size_t i = 0; i < npart; i++)
	{
		double pos = (double(i) - (npart / 2)) * 0.1 * millimeter;
		Particle p1(0);
		p1.x() = pos;
		pcoords.push_back(p1);

		Particle p2(0);
		p2.y() = pos;
		pcoords.push_back(p2);

		Particle p3(0);
		p3.x() = pos;
		p3.y() = pos;
		pcoords.push_back(p3);

		Particle p4(0);
		p4.x() = pos;
		p4.y() = -pos;
		pcoords.push_back(p4);
	}

	const string fname = "basic_hollow_electron_lens_test_out.dat";
	ofstream of(fname);
	if(!of)
	{
		cerr << "Could not open: " << fname << endl;
		exit(1);
	}
	of << "# case id x xp y yp" << endl;

	for(int testcase = 0; testcase < 6; testcase++)
	{
		vector<Particle> bcoords {pcoords};
		ProtonBunch* myBunch = new ProtonBunch(beam_energy, 1, bcoords);

		int nturns = 1;

		HollowElectronLens *hel;

		cout << "Running case " << testcase << endl;
		switch(testcase)
		{
		case (0):
		{
			hel = new HollowElectronLens("hel1", 0, 0, 5, 0.195, 2.334948339E4, 3.0);
			// 1 = opposite to protons (focussing)
			hel->SetElectronDirection(1);
			hel->SetOpMode(DC);
			hel->SetRadii(2 * millimeter, 5 * millimeter);
			break;
		}
		case (1):
		{
			hel = new HollowElectronLens("hel1", 0, 0, 5, 0.5, 1.5E4, 2.0);
			hel->SetElectronDirection(1);
			// radial (measured) or perfect profile
			hel->SetRadialProfile();
			// for LHC hardware we need to scale the radial profile
			hel->SetLHCRadialProfile();
			hel->SetOpMode(DC);
			hel->SetRadii(2 * millimeter, 5 * millimeter);
			break;
		}
		case (2):
		{
			hel = new HollowElectronLens("hel1", 0, 0, 5, 0.195, 2.334948339E4, 3.0);
			hel->SetElectronDirection(0);
			hel->SetOpMode(DC);
			hel->SetRadii(2 * millimeter, 5 * millimeter);
			break;
		}
		case (3):
		{
			hel = new HollowElectronLens("hel1", 0, 0, 5, 0.7, 2.334948339E4, 3.0);
			hel->SetElectronDirection(1);
			hel->SetAC(0.31, .002, 5E-5, 1E3, 2.);
			hel->SetOpMode(AC);
			hel->SetRadii(3 * millimeter, 5 * millimeter);
			break;
		}
		case (4):
		{
			hel = new HollowElectronLens("hel1", 0, 0, 5, 0.195, 2.0E4, 3.0);
			hel->SetElectronDirection(0);
			hel->SetAC(0.31, .002, 5E-5, 1E3, 2.);
			hel->SetOpMode(AC);
			hel->SetRadii(3 * millimeter, 5 * millimeter);
			nturns = 5;
			break;
		}
		case (5):
		{
			hel = new HollowElectronLens("hel1", 0, 0, 5, 0.195, 2.334948339E4, 2.0);
			hel->SetTurnskip(2);
			hel->SetOpMode(Turnskip);
			hel->SetRadii(1 * millimeter, 4 * millimeter);
			nturns = 5;
			break;
		}
		default:
		{
			cerr << "Unknown test case" << endl;
			abort();
		}
		}

		AcceleratorModelConstructor* ctor = new AcceleratorModelConstructor();
		ctor->NewModel();
		ctor->AppendComponent(*hel);
		AcceleratorModel* theModel = ctor->GetModel();
		delete ctor;
		AcceleratorModel::RingIterator ring = theModel->GetRing();
		ParticleTracker* tracker = new ParticleTracker(ring, myBunch);

		auto nancheck = new NANCheckProcess;
		nancheck->SetDetailed(0);
		nancheck->SetHaltNAN(1);
		tracker->AddProcess(nancheck);

		HollowELensProcess* myHELProcess = new HollowELensProcess(3);
		tracker->AddProcess(myHELProcess);

		for(int n = 0; n < nturns; n++)
		{
			cout << "turn: " << n << endl;
			tracker->Track(myBunch);
		}
		size_t i = 0;
		for(auto ip = myBunch->begin(); ip != myBunch->end(); ip++, i++)
		{
			of << testcase << " " << i << " " <<   ip->x() << " " << ip->xp() << " " <<   ip->y() << " " << ip->yp()
			   << endl;
		}

		delete tracker;
		delete myBunch;
		delete theModel;
	}

	of.close();
	return 0;
}
