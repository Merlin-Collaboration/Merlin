/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef utils_h
#define utils_h 1

#include "merlin_config.h"
#include <ctime>
#include <limits>
#include <cmath>

inline bool fequal(double x, double y, double tol = std::numeric_limits<double>::epsilon())
{
	return fabs(x - y) <= tol;
}

inline int Round(double x)
{
	return static_cast<int>(x + 0.5);
}

/**
 * TIMING macro. Used to output the real time used (in seconds)
 * by a function call. The result is output to OS, which must
 * be an ostream.
 */
#define TIMING(FUNC, OS) \
	{time_t start_t = time(0); FUNC; \
	 OS << "done: real time: " << int(difftime(time(0), start_t)) << " seconds" << endl;}

// special function definitions
// ----------------------------

double NormalBin(double x1, double x2);

// Modified Bessel function
double BesselI0(double x);
double BesselI1(double x);
double BesselIn(int n, double x);

// Gamma function
double LogGamma(double);
inline double Gamma(double x)
{
	return exp(LogGamma(x));
}

#endif
