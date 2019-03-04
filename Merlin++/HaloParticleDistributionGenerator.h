/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef HaloParticleDistributionGenerator_h
#define HaloParticleDistributionGenerator_h 1

#include "PSvector.h"
#include "ParticleDistributionGenerator.h"

class Halo1ParticleDistributionGenerator: public ParticleDistributionGenerator
{
public:
	/// @parameter halo_size_ Size of halo in units of sigma
	Halo1ParticleDistributionGenerator(double halo_size_ = 1.0) :
		halo_size(halo_size_)
	{
	}
protected:
	double halo_size;
};

/**
 * Generator for a horizontal halo distribution.
 *
 * Particles in a ring in x, x' at amplitude of halo_size_ sigma. Particles on 0,0 in y, y'
 */
class HorizonalHalo1ParticleDistributionGenerator: public Halo1ParticleDistributionGenerator
{
public:
	using Halo1ParticleDistributionGenerator::Halo1ParticleDistributionGenerator;
	virtual PSvector GenerateFromDistribution() const override;
};

/**
 * Generator for a vertical halo distribution.
 *
 * Particles in a ring in y, y' at amplitude of halo_size_ sigma. Particles on 0,0 in x, x'
 */
class VerticalHalo1ParticleDistributionGenerator: public Halo1ParticleDistributionGenerator
{
public:
	using Halo1ParticleDistributionGenerator::Halo1ParticleDistributionGenerator;
	virtual PSvector GenerateFromDistribution() const override;
};

class Halo2ParticleDistributionGenerator: public ParticleDistributionGenerator
{
public:
	Halo2ParticleDistributionGenerator(double halo_size_ = 1.0, PSvector cutoffs_ = PSvector(0)) :
		halo_size(halo_size_), cutoffs(cutoffs_)
	{
	}
	Halo2ParticleDistributionGenerator(PSvector cutoffs_) :
		halo_size(1.0), cutoffs(cutoffs_)
	{
	}
protected:
	double halo_size;
	PSvector cutoffs;
};

/**
 * Generator for a horizontal halo distribution.
 *
 * Particles in a ring in x, x' at amplitude of halo_size_ sigma. Particles normally distributed y, y'
 */
class HorizonalHalo2ParticleDistributionGenerator: public Halo2ParticleDistributionGenerator
{
public:
	using Halo2ParticleDistributionGenerator::Halo2ParticleDistributionGenerator;
	virtual PSvector GenerateFromDistribution() const override;
};

/**
 * Generator for a vertical halo distribution.
 *
 * Particles in a ring in y, y' at amplitude of halo_size_ sigma. Particles normally distributed in x, x'
 */
class VerticalHalo2ParticleDistributionGenerator: public Halo2ParticleDistributionGenerator
{
public:
	using Halo2ParticleDistributionGenerator::Halo2ParticleDistributionGenerator;
	virtual PSvector GenerateFromDistribution() const override;
};

#endif
