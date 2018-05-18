/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _Normal_h
#define _Normal_h

#include "Random.h"

/**
*
*  See Simulation, Modelling & Analysis by Law & Kelton, pp259
*
*  This is the ``polar'' method.
*
*/

class Normal: public Random
{
	char haveCachedNormal;
	double cachedNormal;

protected:
	double pMean;
	double pVariance;
	double pStdDev;

public:
	Normal(double xmean, double xvariance, RNG *gen);
	double mean();
	double mean(double x);
	double variance();
	double variance(double x);
	virtual double operator()();
};


inline Normal::Normal(double xmean, double xvariance, RNG *gen)
	: Random(gen)
{
	pMean = xmean;
	pVariance = xvariance;
	pStdDev = sqrt(pVariance);
	haveCachedNormal = 0;
}

inline double Normal::mean()
{
	return pMean;
}
inline double Normal::mean(double x)
{
	double t=pMean;
	pMean = x;
	return t;
}

inline double Normal::variance()
{
	return pVariance;
}
inline double Normal::variance(double x)
{
	double t=pVariance;
	pVariance = x;
	pStdDev = sqrt(pVariance);
	return t;
}

#endif
