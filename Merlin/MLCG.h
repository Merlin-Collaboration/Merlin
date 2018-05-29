/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _MLCG_h
#define _MLCG_h 1

#include "RNG.h"
#include <cmath>

/**
*
*	Multiplicative Linear Congruential Generator
*
*/

class MLCG : public RNG
{
	_G_int32_t initialSeedOne;
	_G_int32_t initialSeedTwo;
	_G_int32_t seedOne;
	_G_int32_t seedTwo;

protected:

public:
	MLCG(_G_int32_t seed1 = 0, _G_int32_t seed2 = 1);
	/**
	*
	* Return a long-words word of random bits
	*
	*/
	virtual _G_uint32_t asLong();
	virtual void reset();
	_G_int32_t seed1();
	void seed1(_G_int32_t);
	_G_int32_t seed2();
	void seed2(_G_int32_t);
	void reseed(_G_int32_t, _G_int32_t);
};

inline _G_int32_t
MLCG::seed1()
{
	return seedOne;
}

inline void
MLCG::seed1(_G_int32_t s)
{
	initialSeedOne = s;
	reset();
}

inline _G_int32_t
MLCG::seed2()
{
	return seedTwo;
}

inline void
MLCG::seed2(_G_int32_t s)
{
	initialSeedTwo = s;
	reset();
}

inline void
MLCG::reseed(_G_int32_t s1, _G_int32_t s2)
{
	initialSeedOne = s1;
	initialSeedTwo = s2;
	reset();
}

#endif
