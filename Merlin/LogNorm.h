/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _LogNormal_h
#define _LogNormal_h

#include "Normal.h"

class LogNormal: public Normal
{
protected:
	double logMean;
	double logVariance;
	void setState();
public:
	LogNormal(double mean, double variance, RNG *gen);
	double mean();
	double mean(double x);
	double variance();
	double variance(double x);
	virtual double operator()();
};

inline void LogNormal::setState()
{
	double m2 = logMean * logMean;
	pMean = log(m2 / sqrt(logVariance + m2));
	// from ch@heike.informatik.uni-dortmund.de:
	// (was   pVariance = log((sqrt(logVariance + m2)/m2 )); )
	pStdDev = sqrt(log((logVariance + m2) / m2));
}

inline LogNormal::LogNormal(double mean, double variance, RNG *gen) :
	Normal(mean, variance, gen)
{
	logMean = mean;
	logVariance = variance;
	setState();
}

inline double LogNormal::mean()
{
	return logMean;
}

inline double LogNormal::mean(double x)
{
	double t = logMean;
	logMean = x;
	setState();
	return t;
}

inline double LogNormal::variance()
{
	return logVariance;
}

inline double LogNormal::variance(double x)
{
	double t = logVariance;
	logVariance = x;
	setState();
	return t;
}

#endif
