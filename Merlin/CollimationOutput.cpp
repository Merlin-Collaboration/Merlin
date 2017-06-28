/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 5.01 (2015)
//
// Copyright: see Merlin/copyright.txt
//
// Created:		08.01.15 Haroon Rafique
// Modified:	09.11.15 Haroon Rafique & Alessia Valloni
// Last Edited: 10.11.15 HR + AV
//
/////////////////////////////////////////////////////////////////////////
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

