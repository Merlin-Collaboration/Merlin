/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef Dispersion_h
#define Dispersion_h 1

#include <fstream>
#include <iostream>
#include "AcceleratorModel.h"

class Dispersion
{
public:
	Dispersion(AcceleratorModel* aModel, double refMomentum);
	void FindDispersion(int n = 0);
	void FindRMSDispersion(ofstream* file = nullptr);
	double Dx, Dxp, Dy, Dyp;
	double DxRMS, DyRMS;
	double SetDelta(double new_delta);

private:
	AcceleratorModel* theModel;
	double p0;
	double delta;
};

#endif
