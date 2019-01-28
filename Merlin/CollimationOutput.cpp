/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>

#include "CollimationOutput.h"
#include "MerlinException.h"

namespace ParticleTracking
{

CollimationOutput::CollimationOutput(OutputType ot)
{
	otype = ot;
}

} // End namespace ParticleTracking
