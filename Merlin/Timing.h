/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_Timing
#define _h_Timing 1

#include <ctime>
#include <iostream>

#define BEGIN_TIMING { time_t start_t = time(0)
#define END_TIMING(os,txt) os<<txt<<difftime(time(0),start_t)<<" seconds"; }

#endif
