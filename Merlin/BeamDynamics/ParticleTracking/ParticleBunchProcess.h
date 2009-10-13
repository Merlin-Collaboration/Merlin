/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef ParticleBunchProcess_h
#define ParticleBunchProcess_h 1

#include "merlin_config.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include "BeamDynamics/BunchProcess.h"

namespace ParticleTracking {

typedef TBunchProc<ParticleBunch> ParticleBunchProcess;

};

#endif
