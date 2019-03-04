/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ConstructSrot
#define _h_ConstructSrot 1

#include <string>
#include "ComponentFrame.h"

ComponentFrame* ConstructSrot(double angle, const std::string& name);
ComponentFrame* ConstructXrot(double angle, const std::string& name);

#endif
