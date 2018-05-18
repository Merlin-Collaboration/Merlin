/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _DiscreteUniform_h
#define _DiscreteUniform_h 1

#include "Random.h"

//
//	The interval [lo..hi)
//

class DiscreteUniform: public Random
{
	long pLow;
	long pHigh;
	double delta;
public:
	DiscreteUniform(long low, long high, RNG *gen);

	long low();
	long low(long x);
	long high();
	long high(long x);

	virtual double operator()();
};


inline DiscreteUniform::DiscreteUniform(long low, long high, RNG *gen)
	: Random(gen)
{
	pLow = (low < high) ? low : high;
	pHigh = (low < high) ? high : low;
	delta = (pHigh - pLow) + 1;
}

inline long DiscreteUniform::low()
{
	return pLow;
}

inline long DiscreteUniform::low(long x)
{
	long tmp = pLow;
	pLow = x;
	delta = (pHigh - pLow) + 1;
	return tmp;
}

inline long DiscreteUniform::high()
{
	return pHigh;
}

inline long DiscreteUniform::high(long x)
{
	long tmp = pHigh;
	pHigh = x;
	delta = (pHigh - pLow) + 1;
	return tmp;
}

#endif
