/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 * This file is derived from software bearing the copyright notice: Copyright (c) 1999 by N.J.Walker.  ALL RIGHTS RESERVED.
 */

#ifndef ParticleMap_h
#define ParticleMap_h 1

#include "merlin_config.h"

namespace ParticleTracking
{

class ParticleBunch;

/**
 *	An arbitrary map that can be applied to a ParticleBunch.
 */
class ParticleMap
{
public:

	virtual ~ParticleMap();

	/**
	 *	Applies the map to the specified ParticleBunch.
	 */
	virtual ParticleBunch& Apply(ParticleBunch& bunch) const = 0;

	virtual void Invert() = 0;

protected:
private:
private:
};

/**
 * Class ParticleMap
 */
inline ParticleMap::~ParticleMap()
{

}

} //end namespace ParticleTracking

#endif
