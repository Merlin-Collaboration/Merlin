/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _BUNCHFILTER_H_
#define _BUNCHFILTER_H_

#include "PSvector.h"

namespace ParticleTracking
{

class ParticleBunchFilter
{
public:

	virtual ~ParticleBunchFilter();

	/**
	 *	Used by the ParticleBunch constructor object to select
	 *	vectors for inclusion in a ParticleBunch.
	 */
	virtual bool Apply(const PSvector& v) const = 0;
};

class HorizontalHaloParticleBunchFilter: public ParticleBunchFilter
{
public:

	/**
	 *	Used by the ParticleBunch constructor object to select
	 *	vectors for inclusion in a ParticleBunch.
	 */
	bool Apply(const PSvector& v) const;

	void SetHorizontalLimit(double);
	void SetHorizontalOrbit(double);

private:
	double limit;
	double orbit;
};

class VerticalHaloParticleBunchFilter: public ParticleBunchFilter
{
public:

	/**
	 *	Used by the ParticleBunch constructor object to select
	 *	vectors for inclusion in a ParticleBunch.
	 */
	bool Apply(const PSvector& v) const;

	void SetVerticalLimit(double);

private:
	double limit;
};

}   //End particle tracking namespace

#endif
