/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef LinearParticleMap_h
#define LinearParticleMap_h 1

#include "merlin_config.h"
#include "MatrixMaps.h"
#include "ParticleMap.h"

namespace ParticleTracking
{

class ParticleBunch;

class LinearParticleMap: public ParticleMap
{
public:

	/**
	 *	Applies the map to the specified ParticleBunch.
	 */
	virtual ParticleBunch& Apply(ParticleBunch& bunch) const;
	virtual void Invert();
	RMtrx R;
};

} //end namespace ParticleTracking

#endif
