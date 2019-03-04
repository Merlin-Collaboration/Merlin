/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <sstream>
#include "../tests.h"
#include "ParticleBunchTypes.h"

using namespace std;

int main(int argc, char* argv[])
{
	double beam_mom = 100;
	double charge = 1;

	ParticleBunch* b1 = new ParticleBunch(beam_mom, charge);
	ParticleBunch* b2 = new ParticleBunch(beam_mom, charge);
	assert(b1->size() == 0);
	assert(b2->size() == 0);

	const int npart = 100;
	for(int i = 0; i < npart; i++)
	{
		Particle p(0);
		p.x() = i * 0.1;
		p.xp() = i * 0.1 + 0.01;
		p.y() = i * 0.2;
		p.yp() = i * -0.1;

		b1->AddParticle(p);
	}

	assert(b1->size() == npart);
	assert(b2->size() == 0);

	stringstream bs;

	b1->Output(bs);
	//cout << bs.str() << endl;
	b2->Input(1, bs);

	assert(b1->size() == npart);
	assert(b2->size() == npart);

	for(int i = 0; i < npart; i++)
	{
		Particle &p1 = b1->GetParticles()[i];
		Particle &p2 = b2->GetParticles()[i];
		assert_close(p1.x(), p2.x(), 1e-20);
		assert_close(p1.xp(), p2.xp(), 1e-20);
		assert_close(p1.y(), p2.y(), 1e-20);
		assert_close(p1.yp(), p2.yp(), 1e-20);
	}

	delete b1;
	delete b2;
	return 0;
}
