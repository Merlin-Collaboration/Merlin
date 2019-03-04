/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef CollimatorPotentialModels_h
#define CollimatorPotentialModels_h 1

#include "merlin_config.h"
#include "CollimatorWakePotentials.h"

/**
 * The geometric wake potential:
 * Steeply tapered collimator moving from aperture b to aperture a
 * Ref.:
 * B.W.Zotter and S.A.Kheifets, Impedances and Wakes in High-Energy Particle Accelerators,
 * World Scientific (1998)
 */
class TaperedCollimatorPotentials: public CollimatorWakePotentials
{
public:
	TaperedCollimatorPotentials(int m, double aa, double bb);
	~TaperedCollimatorPotentials();
	virtual double Wlong(double z, int m) const;
	virtual double Wtrans(double z, int m) const;
	double Wlong(double z) const
	{
		return 0;
	}
	double Wtrans(double z) const
	{
		return 0;
	}
private:
	double* coeff;
	double a, b;
};

/**
 * the resistive wake potentials  (in MKS system)
 */
class ResistiveWakePotentials: public CollimatorWakePotentials
{
public:
	ResistiveWakePotentials(int m, double r, double s, double l);
	~ResistiveWakePotentials();
	virtual double Wlong(double z, int m) const;
	virtual double Wtrans(double z, int m) const;

private:
	double* coeff;
	double rad, sigma, length;
};
#endif
