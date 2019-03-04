/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CollimatorWakePotentials_h
#define CollimatorWakePotentials_h 1

#include "merlin_config.h"

#include "WakePotentials.h"

/**
 * Abstract class for calculating the longitudinal and
 * transverse single-bunch wakefield potentials (Greens
 * functions) with modes
 */
class CollimatorWakePotentials: public WakePotentials
{

public:

	CollimatorWakePotentials(int m, double rad = 0, double conduct = 0)
	//take the radius and the conductivity out of WakePotentials
		: WakePotentials(rad, conduct), nmodes(m)
	{
	}

	virtual ~CollimatorWakePotentials()
	{
	}

	virtual double Wlong(double s, int m) const = 0;
	virtual double Wtrans(double s, int m) const = 0;

protected:

	int nmodes;

private:

	using WakePotentials::Wlong;
	using WakePotentials::Wtrans;

};

#endif
