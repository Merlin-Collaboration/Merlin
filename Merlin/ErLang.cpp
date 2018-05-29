/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#include "Random.h"
#include "Erlang.h"

double Erlang::operator()()
{
	double prod = 1.0;

	for (int i = 0; i < k; i++)
	{
		prod *= pGenerator -> asDouble();
	}
	return -log(prod)/a;
}
