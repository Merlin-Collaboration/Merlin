/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#include "Random.h"
#include "Binomial.h"

double Binomial::operator()()
{
	int s = 0;
	for(int i = 0; i < pN; i++)
	{
		if(pGenerator->asDouble() < pU)
		{
			s++;
		}
	}
	return double(s);
}
