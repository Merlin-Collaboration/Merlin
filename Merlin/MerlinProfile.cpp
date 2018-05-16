/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "MerlinProfile.h"

#include <map>
#include <string>
std::map<std::string, timespec> MerlinProfile::pData;
std::map<std::string, timespec> MerlinProfile::StartTime;

// To test whether MERLIN_PROFILE was enabled when libmerlin was built
#ifdef MERLIN_PROFILE
bool MerlinProfile::IsEnabled()
{
	return true;
}
#else
bool MerlinProfile::IsEnabled()
{
	return false;
}
#endif
