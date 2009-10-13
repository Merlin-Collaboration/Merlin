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

#ifndef ParticleMapPI_h
#define ParticleMapPI_h 1

#include "merlin_config.h"
// ParticleBunchIntegrator
#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"
// ParticleMapComponent
#include "BeamDynamics/ParticleTracking/ParticleMapComponent.h"

namespace ParticleTracking {

class ParticleMapCI : public ParticleComponentTracker::Integrator<ParticleMapComponent> {
public:
    virtual void TrackStep (double ds);
};

}; // end namespace ParticleTracking

#endif
