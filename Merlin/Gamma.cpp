/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "utils.h"

double LogGamma(double xx)
{
	double x, tmp, ser;
	static double cof[6] = {76.18009173, -86.50532033, 24.01409822,
							-1.231739516, 0.120858003e-2, -0.536382e-5};
	int j;

	x = xx - 1.0;

	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.0;
	for(j = 0; j <= 5; j++)
	{
		x += 1.0;
		ser += cof[j] / x;
	}
	return -tmp + log(2.50662827465 * ser);
}
