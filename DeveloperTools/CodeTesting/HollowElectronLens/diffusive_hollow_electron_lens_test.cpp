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
#include "HollowElectronLens.h"
#include "HollowELensProcess.h"
#include "RandomNG.h"
#include "NANCheckProcess.h"

/*
 * Test that the diffusive modes of the the Hollow Election Lens matches
 * the DC mode after an equivalent number of passes. I.E checks that the
 * diffusive mode fires on approximately 50% of passes (2 sigma check).
 * Note that this is just a single HEL, so no phase advance between turns.
 *
 */

using namespace std;
using namespace PhysicalUnits;

int main()
{
	HollowElectronLens *diff_hel = new HollowElectronLens("hel1", 0, 2, 5, 0.195, 2.334948339E4, 3.0);
	diff_hel->SetElectronDirection(1);
	diff_hel->SetRadii(2 * millimeter, 5 * millimeter);

	HollowElectronLens *ac_hel = new HollowElectronLens("hel1", 0, 0, 5, 0.195, 2.334948339E4, 3.0);
	ac_hel->SetElectronDirection(1);
	ac_hel->SetRadii(2 * millimeter, 5 * millimeter);

	const double beam_energy = 7000.0;
	const size_t npart = 400;
	const int diff_turns = 100;

	int seed = 1; // seed 0 starts with run of 11 negative numbers
	RandomNG::init(seed);

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

	// first track diff_turns turns though the diffusive mode

	vector<Particle> diff_coords {pcoords};
	ProtonBunch* diff_bunch = new ProtonBunch(beam_energy, 1, diff_coords);

	AcceleratorModelConstructor* ctor = new AcceleratorModelConstructor();
	ctor->NewModel();
	ctor->AppendComponent(*diff_hel);
	AcceleratorModel* theModel = ctor->GetModel();
	delete ctor;

	AcceleratorModel::RingIterator ring = theModel->GetRing();
	ParticleTracker* diff_tracker = new ParticleTracker(ring, diff_bunch);

	auto nancheck = new NANCheckProcess;
	nancheck->SetDetailed(0);
	nancheck->SetHaltNAN(1);
	diff_tracker->AddProcess(nancheck);

	HollowELensProcess* HELProcess = new HollowELensProcess(3);
	diff_tracker->AddProcess(HELProcess);

	for(int n = 0; n < diff_turns; n++)
	{
		diff_tracker->Track(diff_bunch);
	}
	//diff_bunch->begin()->y() += 1e-6; // artificial error to test the test
	delete theModel;
	delete diff_tracker;

	// then track though the DC mode, until the bunch matches the diffusive
	vector<Particle> ac_coords {pcoords};
	ProtonBunch* ac_bunch = new ProtonBunch(beam_energy, 1, ac_coords);

	ctor = new AcceleratorModelConstructor();
	ctor->NewModel();
	ctor->AppendComponent(*ac_hel);
	theModel = ctor->GetModel();
	delete ctor;
	ring = theModel->GetRing();
	ParticleTracker* ac_tracker = new ParticleTracker(ring, ac_bunch);

	auto nancheck2 = new NANCheckProcess;
	nancheck2->SetDetailed(0);
	nancheck2->SetHaltNAN(1);
	ac_tracker->AddProcess(nancheck2);

	HELProcess = new HollowELensProcess(3);
	ac_tracker->AddProcess(HELProcess);

	int found_match = -1;
	for(int n = 0; n < diff_turns; n++)
	{
		ac_tracker->Track(ac_bunch);

		bool match = true;
		size_t i = 0;
		for(auto ipd = diff_bunch->begin(), ipa = ac_bunch->begin(); ipd != diff_bunch->end(); ipd++, ipa++, i++)
		{
			//cout << " " << i << " " <<   ipd->x() << " " << ipd->xp() << " " <<   ipd->y() << " " << ipd->yp() <<endl;
			//cout << " " << i << " " <<   ipa->x() << " " << ipa->xp() << " " <<   ipa->y() << " " << ipa->yp() <<endl;

			if(ipd->x() != ipa->x()
				|| ipd->xp() != ipa->xp()
				|| ipd->y() != ipa->y()
				|| ipd->yp() != ipa->yp()
				)
			{
				match = false;
			}
		}
		if(match)
		{
			found_match = n;
			break;
		}
	}

	delete ac_tracker;
	delete diff_bunch;
	delete ac_bunch;
	delete theModel;

	if(found_match != -1)
	{
		cout << "Found match on turn " << found_match << endl;
		cout << "(Expected  " << diff_turns / 2.0 << ")" << endl;
	}
	else
	{
		cout << "No match found" << endl;
		return 1;
	}
	// check within 2 sigma
	if(fabs(found_match - diff_turns / 2.0) < 2 * sqrt(diff_turns / 2.0))
	{
		cout << "Pass" << endl;
		return 0;
	}
	else
	{
		cout << "Fail" << endl;
		return 1;
	}
}
