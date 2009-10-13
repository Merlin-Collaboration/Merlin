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

#include "BeamDynamics/ParticleTracking/Integrators/ParticleMapPI.h"

namespace ParticleTracking {

void ParticleMapCI::TrackStep (double ds)
{
    assert(ds==0);
    currentComponent->Apply(*currentBunch);
    return;
}

}; // end namespace ParticleTracking
