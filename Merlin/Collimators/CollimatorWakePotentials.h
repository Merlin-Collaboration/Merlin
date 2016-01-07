/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 3 (2004)
//
// Copyright: see Merlin/copyright.txt
//
//  Created: June, 2006
//
/////////////////////////////////////////////////////////////////////////

#ifndef CollimatorWakePotentials_h
#define CollimatorWakePotentials_h 1

#include "merlin_config.h"

#include "AcceleratorModel/WakePotentials.h"

//	Abstract class for calculating the longitudinal and
//	transverse single-bunch wakefield potentials (Greens
//	functions) with modes

class CollimatorWakePotentials : public WakePotentials
{


public:

	CollimatorWakePotentials(int m, double r=0, double s=0)
		: WakePotentials(r,s)
	{
		nmodes = m;    //take the radius and the conductivity out of WakePotentials
	}

	virtual ~CollimatorWakePotentials () {};

	virtual double Wlong (double s, int m) const = 0;
	virtual double Wtrans(double s, int m) const = 0;

protected:

	int nmodes;
};

#endif
