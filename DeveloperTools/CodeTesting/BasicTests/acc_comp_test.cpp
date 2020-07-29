/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>
#include <memory>

#include "MerlinVersion.h"
#include "Components.h"

/*
 * Simple tests creating elements, getting and setting parameters
 */

using namespace std;

void basic_compents()
{
	vector<unique_ptr<const AcceleratorComponent> > comps;
	comps.push_back(make_unique<Drift>("d", 1));
	comps.push_back(make_unique<SectorBend>("b1", 2, 100, 2));
	comps.push_back(make_unique<Quadrupole>("q1", 3, 20));
	comps.push_back(make_unique<SkewQuadrupole>("sq1", 4, 20));
	comps.push_back(make_unique<Sextupole>("sex1", 5, 20));
	comps.push_back(make_unique<SkewSextupole>("ssex1", 6, 20));
	comps.push_back(make_unique<Octupole>("oct1", 7, 20));
	comps.push_back(make_unique<Decapole>("dec1", 8, 20));
	comps.push_back(make_unique<Marker>("m1"));

	for(const auto &el : comps)
	{
		cout << el->GetType() << " " << el->GetName() << " " << el->GetLength() << endl;
	}
	assert(comps[1]->GetName() == "b1");
	assert(comps[1]->GetQualifiedName() == "SectorBend.b1");
	assert(comps[2]->GetQualifiedName() == "Quadrupole.q1");
}

void modify_quad()
{
	auto q1 = Quadrupole("quad", 3, 10);

	assert(q1.GetLength() == 3);
	assert(q1.GetFieldStrength() == 10);

	assert(q1.GetEMField()->GetBFieldAt({1, 0, 0}, 0).y == 10);
	assert(q1.GetEMField()->GetBFieldAt({0.5, 0, 0}, 0).y == 5);
	assert(q1.GetEMField()->GetBFieldAt({0, 1, 0}, 0).y == 0);
	assert(q1.GetEMField()->GetBFieldAt({1, 0, 0}, 0).x == 0);
	assert(q1.GetEMField()->GetBFieldAt({0, 1, 0}, 0).x == 10);

	q1.SetFieldStrength(12);
	assert(q1.GetFieldStrength() == 12);
	assert(q1.GetEMField()->GetBFieldAt({1, 0, 0}, 0).y == 12);

	q1.SetName("f1");
	assert(q1.GetName() == "f1");

}

int main(int argc, char* argv[])
{
	cout << merlin_version_info();

	basic_compents();
	modify_quad();

	cout << "All element tests successful" << endl;
}
