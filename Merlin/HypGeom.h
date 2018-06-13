/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _HyperGeometric_h
#define _HyperGeometric_h

#include "Random.h"

class HyperGeometric: public Random
{
protected:
	double pMean;
	double pVariance;
	double pP;
	void setState();

public:
	HyperGeometric(double mean, double variance, RNG *gen);

	double mean();
	double mean(double x);
	double variance();
	double variance(double x);

	virtual double operator()();
};

inline void HyperGeometric::setState()
{
	double z = pVariance / (pMean * pMean);
	pP = 0.5 * (1.0 - sqrt((z - 1.0) / (z + 1.0)));
}

inline HyperGeometric::HyperGeometric(double mean, double variance, RNG *gen) :
	Random(gen)
{
	pMean = mean;
	pVariance = variance;
	setState();
}

inline double HyperGeometric::mean()
{
	return pMean;
}

inline double HyperGeometric::mean(double x)
{
	double t = pMean;
	pMean = x;
	setState();
	return t;
}

inline double HyperGeometric::variance()
{
	return pVariance;
}

inline double HyperGeometric::variance(double x)
{
	double t = pVariance;
	pVariance = x;
	setState();
	return t;
}

#endif
