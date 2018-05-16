/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _Uniform_h
#define _Uniform_h 1

#include "Random.h"

/**
*	The interval [lo..hi]
*/
class Uniform: public Random
{
	double pLow;
	double pHigh;
	double delta;
public:
	Uniform(double low, double high, RNG *gen);

	double low();
	double low(double x);
	double high();
	double high(double x);

	virtual double operator()();
};


inline Uniform::Uniform(double low, double high, RNG *gen) : Random(gen)
{
	pLow = (low < high) ? low : high;
	pHigh = (low < high) ? high : low;
	delta = pHigh - pLow;
}

inline double Uniform::low()
{
	return pLow;
}

inline double Uniform::low(double x)
{
	double tmp = pLow;
	pLow = x;
	delta = pHigh - pLow;
	return tmp;
}

inline double Uniform::high()
{
	return pHigh;
}

inline double Uniform::high(double x)
{
	double tmp = pHigh;
	pHigh = x;
	delta = pHigh - pLow;
	return tmp;
}

#endif
