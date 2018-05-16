/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>

#define assert_close(x1, x2, d) assert( fabs(x1 - x2) < d )

#define assert_throws(func_call, exception)\
try{(func_call); fprintf(stderr, "Assertion "#func_call" throws "#exception" failed: Nothing thrown\n"); abort();}\
catch(exception &e){}\
catch(...){fprintf(stderr, "Assertion "#func_call" throws "#exception" failed: Other exception thrown\n");}

