/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "../tests.h"

#include <iostream>

#include "RandomNG.h"
#include "AcceleratorModelConstructor.h"
#include "Components.h"
#include "SupportStructure.h"
#include "ATL2D.h"
#include "SimpleATL.h"

/*
 * Run a few simple cases of the ATL2D and SimpleATL
 */

int main()
{
	RandomNG::init();

	AcceleratorModelConstructor am_ctor;
	am_ctor.NewModel();

	GirderMount* g1 = new GirderMount("g1");
	GirderMount* g2 = new GirderMount("g2");

	Drift* d1 = new Drift("d1", 1);
	Drift* d2 = new Drift("d2", 1);
	Drift* d3 = new Drift("d3", 1);
	Drift* d4 = new Drift("d4", 1);

	Quadrupole* q1 = new Quadrupole("q1", 0.5, 0.2);
	Quadrupole* q2 = new Quadrupole("q1", 0.5, -0.2);

	am_ctor.NewFrame(g1);
	am_ctor.AppendComponent(*d1);
	am_ctor.AppendComponent(*q1);
	am_ctor.AppendComponent(*d2);
	am_ctor.EndFrame();

	am_ctor.NewFrame(g2);
	am_ctor.AppendComponent(*d3);
	am_ctor.AppendComponent(*q2);
	am_ctor.AppendComponent(*d4);
	am_ctor.EndFrame();

	AcceleratorModel* model = am_ctor.GetModel();

	AcceleratorSupportList supports;
	size_t ns = model->GetAcceleratorSupports(supports);

	cout << "Number of supports: " << ns << endl;
	cout << "Initial offsets" << endl;
	for(auto sup : supports)
	{
		auto offset = sup->GetOffset();
		cout << " " << offset.x << " " <<  offset.y << " " <<  offset.z << endl;
	}

	ATL2D atl{1e-3, supports};

	cout << "ATL2D offsets absolute" << endl;
	atl.SetRandomSeed(21);
	atl.SetVibration(1e-3);
	atl.DoStep(0.1);
	atl.RecordOffsets(cout);

	for(auto sup : supports)
	{
		auto offset = sup->GetOffset();
		cout << " " << offset.x << " " <<  offset.y << " " <<  offset.z << endl;
		assert(offset.x == 0);
		assert(offset.y != 0);
		assert(offset.z == 0);
	}

	atl.Reset();
	cout << "ATL2D offsets reset" << endl;
	for(auto sup : supports)
	{
		auto offset = sup->GetOffset();
		cout << " " << offset.x << " " <<  offset.y << " " <<  offset.z << endl;
		assert(offset.x == 0);
		assert(offset.y == 0);
		assert(offset.z == 0);
	}

	cout << "ATL2D offsets increment" << endl;
	atl.SetATLMode(ATL2D::increment);
	atl.DoStep(0.1);
	atl.RecordOffsets(cout);
	for(auto sup : supports)
	{
		auto offset = sup->GetOffset();
		cout << " " << offset.x << " " <<  offset.y << " " <<  offset.z << endl;
		assert(offset.x == 0);
		assert(offset.z == 0);
	}
	assert(supports[0]->GetOffset().y == 0);
	assert(supports[1]->GetOffset().y != 0);
	assert(supports[2]->GetOffset().y != 0);
	assert(supports[3]->GetOffset().y != 0);

	SimpleATL satl{1e-3, supports};
	satl.SetRandomSeed(21);
	satl.Reset();
	cout << "SimpleATL offsets reset" << endl;
	for(auto sup : supports)
	{
		auto offset = sup->GetOffset();
		cout << " " << offset.x << " " <<  offset.y << " " <<  offset.z << endl;
		assert(offset.x == 0);
		assert(offset.y == 0);
		assert(offset.z == 0);
	}

	cout << "SimpleATL offsets" << endl;
	satl.DoStep(0.1);
	satl.RecordOffsets(cout);

	for(auto sup : supports)
	{
		auto offset = sup->GetOffset();
		cout << " " << offset.x << " " <<  offset.y << " " <<  offset.z << endl;
		assert(offset.x == 0);
		assert(offset.z == 0);
	}
	assert(supports[0]->GetOffset().y == 0);
	assert(supports[1]->GetOffset().y != 0);
	assert(supports[2]->GetOffset().y != 0);
	assert(supports[3]->GetOffset().y != 0);

	delete model;
	return 0;
}
