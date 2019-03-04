/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_CommonDataStructures
#define _h_CommonDataStructures 1

#include "Klystron.h"
#include "ReferenceParticle.h"
#include <iostream>
#include <vector>
#include <utility>

typedef std::vector<Klystron*> KlystronArray;
typedef std::vector<ReferenceParticle*> ReferenceParticleArray;
typedef std::pair<size_t, size_t> DFS_Segment;
typedef std::vector<size_t> IntegerArray;

inline std::ostream& operator<<(std::ostream& os, const DFS_Segment& seg)
{
	return os << '[' << seg.first << ',' << seg.second << ']';
}

#endif
