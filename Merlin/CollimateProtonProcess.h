/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//
// Class library version 5.01 (2015)
//
// Copyright: see Merlin/copyright.txt
//
// Created:		2010	 RJB
// Modified:	07.09.15 Haroon Rafique
// Last Edited: 03.11.15 HR
//
/////////////////////////////////////////////////////////////////////////
#ifndef CollimateProtonProcess_h
#define CollimateProtonProcess_h 1

#include <fstream>

#include "CollimateParticleProcess.h"
#include "ScatteringModel.h"

using namespace Collimation;

namespace ParticleTracking
{

class CollimateProtonProcess : public CollimateParticleProcess
{
public:

	/**
	* Constructor taking the collimation mode, and the output
	* stream pointer to which to print the results. mode can
	* be a logical OR combination of the collimation modes. A
	* null pointer for osp (default) suppresses output.
	*/
	CollimateProtonProcess (int priority, int mode, std::ostream* osp = nullptr);

	void SetScatteringModel(Collimation::ScatteringModel* s);

private:
	Collimation::ScatteringModel* scattermodel;

	bool DoScatter(Particle&);

};

} // end namespace ParticleTracking

#endif
