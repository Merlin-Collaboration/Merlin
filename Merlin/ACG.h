/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice: Copyright (C) 1988 Free Software Foundation, written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 */

#ifndef _ACG_h
#define _ACG_h 1

#include "RNG.h"
#include <cmath>

/**
 *   \class ACG
 *
 *	Additive number generator. This method is presented in Volume II
 *	of The Art of Computer Programming by Knuth. I've coded the algorithm
 *	and have added the extensions by Andres Nowatzyk of CMU to randomize
 *	the result of algorithm M a bit	by using an LCG & a spatial
 *	permutation table.
 *
 *	The version presented uses the same constants for the LCG that Andres
 *	uses (chosen by trial & error). The spatial permutation table is
 *	the same size (it's based on word size). This is for 32-bit words.
 *
 *	The ``auxillary table'' used by the LCG table varies in size, and
 *	is chosen to be the smallest power of two which is larger than
 *	twice the size of the state table.
 *
 */

typedef unsigned int _G_uint32_t;

class ACG: public RNG
{

	/**
	 * Used to reset generator
	 */
	_G_uint32_t initialSeed;
	int initialTableEntry;

	_G_uint32_t *state;
	_G_uint32_t *auxState;
	short stateSize;
	short auxSize;
	_G_uint32_t lcgRecurr;
	short j;
	short k;

protected:

public:
	ACG(_G_uint32_t seed = 0, int size = 55);
	virtual ~ACG();

	/**
	 * Return a long-words word of random bits
	 * @return A long-words type of random bits
	 */
	virtual unsigned int asLong();
	virtual void reset();

private:
	//Copy protection
	ACG(const ACG& rhs);
	ACG& operator=(const ACG& rhs);
};

#endif
