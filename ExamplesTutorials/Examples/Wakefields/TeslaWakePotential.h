/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _H_TeslaWakePotential
#define _H_TeslaWakePotential

#include "WakePotentials.h"

// A form of WakeFieldProcess::WakePotential that uses linear interpolation
// from a tabulated wakefield file.
class TeslaWakePotentials: public WakePotentials
{
public:

	double Wlong(double z) const;
	double Wtrans(double z) const;
};

#endif
