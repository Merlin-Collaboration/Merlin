//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (1999)
* 
* file Merlin\BeamDynamics\ParticleTracking\ParticleMapComponent.cpp
* last modified 26/09/02 15:18:14
*/

/*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
* Copyright (c) 2001 by The Merlin Collaboration.
* - ALL RIGHTS RESERVED - 
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all copies and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
*/



// ParticleBunch
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
// ComponentTracker
#include "AcceleratorModel/TrackingInterface/ComponentTracker.h"
// ParticleMapComponent
#include "BeamDynamics/ParticleTracking/ParticleMapComponent.h"


#include <cassert>

namespace ParticleTracking {

// Class ParticleMapComponent


const int ParticleMapComponent::ID = UniqueIndex();



ParticleMapComponent::ParticleMapComponent (const std::string& id, ParticleMap* pmap, double intB2ds)

        : AcceleratorComponent(id),itsMap(pmap),ib2(intB2ds)

{

    assert(pmap);

}





const string& ParticleMapComponent::GetType () const
{

    _TYPESTR(ParticleMap);

}


ModelElement* ParticleMapComponent::Copy () const
{

    return new ParticleMapComponent(*this);

}


int ParticleMapComponent::GetIndex () const
{

    return ID;

}


void ParticleMapComponent::RotateY180 ()
{

    itsMap->Invert();

}


void ParticleMapComponent::PrepareTracker (ComponentTracker& aTracker)
{

    _PREPTRACK(aTracker,AcceleratorComponent)

}

}; //end namespace ParticleTracking

