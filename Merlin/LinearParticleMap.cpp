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

#include <cassert>
#include "ParticleBunch.h"
#include "LinearParticleMap.h"

namespace ParticleTracking
{

ParticleBunch& LinearParticleMap::Apply (ParticleBunch& bunch) const
{
	R.Apply(bunch.GetParticles());
	return bunch;
}

void LinearParticleMap::Invert ()
{
	// TODO:
	assert(false);
}

} //end namespace ParticleTracking

