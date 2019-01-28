/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_Histogram
#define _h_Histogram 1

#include <vector>
#include <cstddef>

size_t Hist(const std::vector<double>& data, double x1, double x2, double dx, std::vector<double>& hist);

#endif
