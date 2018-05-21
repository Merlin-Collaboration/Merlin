/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _Random_h
#define _Random_h 1

#include <cmath>
#include "RNG.h"

class Random
{
protected:
	RNG *pGenerator;

public:
	virtual ~Random() {}
	Random(RNG *generator);
	virtual double operator()() = 0;

	RNG *generator();
	void generator(RNG *p);

private:
	//Copy protection
	Random(const Random& rhs);
	Random& operator=(const Random& rhs);
};


inline Random::Random(RNG *gen)
{
	pGenerator = gen;
}

inline RNG *Random::generator()
{
	return(pGenerator);
}

inline void Random::generator(RNG *p)
{
	pGenerator = p;
}

#endif
