/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_model_construction
#define _h_model_construction

#include "merlin_config.h"
#include "AcceleratorModel.h"
#include "ComponentFrame.h"
#include "BeamData.h"
#include <utility>

std::pair<AcceleratorModel*, BeamData*> ConstructModel(const string& fname);

#endif
