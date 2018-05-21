/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#include "Random.h"
#include "Poisson.h"

double Poisson::operator()()
{
	double bound = exp(-1.0 * pMean);
	int count = 0;

	for (double product = 1.0;
	        product >= bound;
	        product *= pGenerator -> asDouble())
	{
		count++;
	}
	return count - 1;
}
