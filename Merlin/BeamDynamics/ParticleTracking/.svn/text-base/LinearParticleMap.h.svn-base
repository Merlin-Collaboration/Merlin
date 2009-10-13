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

#ifndef LinearParticleMap_h
#define LinearParticleMap_h 1

#include "merlin_config.h"
// MatrixMaps
#include "BasicTransport/MatrixMaps.h"
// ParticleMap
#include "BeamDynamics/ParticleTracking/ParticleMap.h"

namespace ParticleTracking {

class ParticleBunch;

class LinearParticleMap : public ParticleMap
{
public:

    //	Applies the map to the specified ParticleBunch.
    virtual ParticleBunch& Apply (ParticleBunch& bunch) const;
    virtual void Invert ();
    RMtrx R;
};

}; //end namespace ParticleTracking

#endif
