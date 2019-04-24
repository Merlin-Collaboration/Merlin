/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef HollowElensParticleDistributionGenerator_h
#define HollowElensParticleDistributionGenerator_h 1

#include "PSvector.h"
#include "ParticleDistributionGenerator.h"

class HELHaloParticleDistributionGenerator: public ParticleDistributionGenerator
{
public:
	/// @parameter halo_inner_ Minimum radius of halo in units of sigma
	/// @parameter halo_outer_ Maximum radius of halo in units of sigma 
	Halo1ParticleDistributionGenerator(double halo_inner_ = 1.0, double halo_outer_ = 1.0 ) :
		halo_inner(halo_inner_) {}
		halo_outer(halo_outer_) {}
protected:
	double halo_inner;
	double halo_outer;
};

/**
 * Generator for a thick HEL halo distribution.
 *
 * Particles in an annulus in x,y between 2 radii halo_inner_ and halo_outer_ in units of sigma, particles at 0,0 in x',y'
 */
class ThickHELHaloParticleDistributionGenerator: public HELHaloParticleDistributionGenerator
{
public:
	using HELHaloParticleDistributionGenerator::HELHaloParticleDistributionGenerator;
	virtual PSvector GenerateFromDistribution() const override;
};

/**
 * Generator for a thin HEL halo distribution.
 *
 * Particles in a ring in x,y at radius halo_outer_ in units of sigma, particles at 0,0 in x',y'
 */
class ThinHELHaloParticleDistributionGenerator: public HELHaloParticleDistributionGenerator
{
public:
	using HELHaloParticleDistributionGenerator::HELHaloParticleDistributionGenerator;
	virtual PSvector GenerateFromDistribution() const override;
};
#endif

