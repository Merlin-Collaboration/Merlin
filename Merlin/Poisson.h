/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _Poisson_h
#define _Poisson_h

#include "Random.h"

class Poisson: public Random
{
protected:
	double pMean;
public:
	Poisson(double mean, RNG *gen);

	double mean();
	double mean(double x);

	virtual double operator()();
};


inline Poisson::Poisson(double mean, RNG *gen)
	: Random(gen)
{
	pMean = mean;
}

inline double Poisson::mean()
{
	return pMean;
}
inline double Poisson::mean(double x)
{
	double t = pMean;
	pMean = x;
	return t;
}

#endif
