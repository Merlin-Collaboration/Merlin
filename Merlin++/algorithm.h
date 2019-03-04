/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef stdext_algorithm_h
#define stdext_algorithm_h 1

#include "merlin_config.h"

template<class InIt, class OutIt, class Pred>
OutIt copy_if(InIt first, InIt last, OutIt x, Pred test)
{
	while(first++ != last)
		if(test(*first))
		{
			*x++ = *first;
		}
	return x;
}

#endif
