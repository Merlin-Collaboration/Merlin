/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "utils.h"
#include <cmath>

using namespace std;

double NormalBin(double x1, double x2)
{
	static const double root2 = sqrt(2.0);
	return 0.5 * (erfc(x1 / root2) - erfc(x2 / root2));
}
