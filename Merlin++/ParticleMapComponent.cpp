/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "ParticleBunch.h"
#include "ComponentTracker.h"
#include "ParticleMapComponent.h"

#include <cassert>

namespace ParticleTracking
{

// Class ParticleMapComponent

const int ParticleMapComponent::ID = UniqueIndex();

ParticleMapComponent::ParticleMapComponent(const std::string& id, ParticleMap* pmap, double intB2ds) :
	AcceleratorComponent(id), itsMap(pmap), ib2(intB2ds)
{
	assert(pmap);
}

const string& ParticleMapComponent::GetType() const
{
	_TYPESTR(ParticleMap);
}

ModelElement* ParticleMapComponent::Copy() const
{
	return new ParticleMapComponent(*this);
}

int ParticleMapComponent::GetIndex() const
{
	return ID;
}

void ParticleMapComponent::RotateY180()
{
	itsMap->Invert();
}

void ParticleMapComponent::PrepareTracker(ComponentTracker& aTracker)
{
	_PREPTRACK(aTracker, AcceleratorComponent)
}
} //end namespace ParticleTracking
