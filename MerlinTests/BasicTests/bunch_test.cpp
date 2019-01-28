/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include "../tests.h"
#include "ParticleBunchTypes.h"

using namespace std;

int main(int argc, char* argv[])
{
	double beam_mom = 100;
	double charge = 1;

	ProtonBunch* myBunch_p = new ProtonBunch(beam_mom, charge);
	ElectronBunch* myBunch_e = new ElectronBunch(beam_mom, charge);

	assert_close(myBunch_p->GetParticleMass(), 1.67262e-27, 1e-31);
	assert_close(myBunch_e->GetParticleMass(), 9.10938215e-31, 1e-38);

	assert(myBunch_p->size() == 0);
	Particle p(0);
	p.x() = 1;
	myBunch_p->AddParticle(p);
	assert(myBunch_p->size() == 1);
	assert(myBunch_p->FirstParticle().x() == 1);

	delete myBunch_p;
	delete myBunch_e;
	return 0;
}
