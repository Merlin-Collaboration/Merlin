/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ConstructModel
#define _h_ConstructModel 1

#include "AcceleratorModel.h"
#include <string>

AcceleratorModel* ConstructModel(const std::string& fname, double E0, const std::string& logFile);

#endif
