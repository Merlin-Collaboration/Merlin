/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>

#include "TeslaWakePotential.h"

//using namespace std;

double TeslaWakePotentials::Wlong(double z) const

{

	return 38.1e+12 * (1.165 * exp(-sqrt(z / 3.65e-3)) - 0.165);

}

double TeslaWakePotentials::Wtrans(double z) const

{

//	return (1290.0e+12*sqrt(z)-2600*z);

//	return 1.0e+12*(1290.0*sqrt(z)-2600.0*z);

	double arZ = sqrt(z / 0.92e-03);

	return 1.21e14 * (1.0 - (1.0 + arZ) * exp(-arZ));

}
