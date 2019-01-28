/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ModelConstruction
#define _h_ModelConstruction

#include <utility>
#include <string>

class AcceleratorModel;
class BeamData;

std::pair<AcceleratorModel*, BeamData*> ConstructModel(const std::string& fname, bool curved);

#endif
