/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _CPUFeatures_h_
#define _CPUFeatures_h_ 1

#include <iostream>
#include <string>

#ifdef LIBNUMA
#include <numa.h>
#endif

namespace CPUFeatures
{

void CheckCPUFeatures();
unsigned int GetCPUFeatures1();
unsigned int GetCPUFeatures2();
std::string GetCPUName();

/**
 * NUMA features (Non Uniform Memory Access)
 *
 * @see `man numa` on linux for details
 */
#ifdef LIBNUMA
bool CheckNUMA();
void PrintNUMAInfo();
#endif

} //End namespace

#endif
