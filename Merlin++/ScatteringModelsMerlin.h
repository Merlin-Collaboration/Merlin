/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ScatteringModelsMerlin_h
#define ScatteringModelsMerlin_h 1

#include "ScatteringModel.h"

// probably not needed, but MaterialProperties added

#include "MaterialProperties.h"

namespace Collimation
{
/**
 * Base class for fixed scattering models
 *
 * Classes derived from this should set up their scattering processes in
 * their constructor, and then set is_fixed to true.
 */
class ScatteringModelFixed: public ScatteringModel
{
public:
	ScatteringModelFixed(MaterialProperties* mat);
	virtual ~ScatteringModelFixed();

	/**
	 * For fixed scattering models, this cannot be called from outside
	 * of the constructor.
	 */
	void AddProcess(Collimation::ScatteringProcess* S);

protected:
	bool is_fixed;
};

/**
 * Merlin physics model
 *
 * ScatteringModel preset with Merlin physics model as described in
 * Appleby, R.B., Barlow, R.J., Molson, J.G. et al.
 * Eur. Phys. J. C (2016) 76: 520.
 * https://dx.doi.org/10.1140/epjc/s10052-016-4363-7
 */
class ScatteringModelMerlin: public ScatteringModelFixed
{
public:
	ScatteringModelMerlin(MaterialProperties* mat);
};

/**
 * Sixtrack style physics model
 *
 * ScatteringModel preset with physics based on SixTrack K2 scattering
 */
class ScatteringModelSixTrack: public ScatteringModelFixed
{
public:
	ScatteringModelSixTrack(MaterialProperties* mat);
};

/**
 * Sixtrack style physics model + new  Ionisation
 */
class ScatteringModelSixTrackIoniz: public ScatteringModelFixed
{
public:
	ScatteringModelSixTrackIoniz(MaterialProperties* mat);
};

/**
 * Sixtrack style physics model + new elastic scattering
 */
class ScatteringModelSixTrackElastic: public ScatteringModelFixed
{
public:
	ScatteringModelSixTrackElastic(MaterialProperties* mat);
};

/**
 * Sixtrack style physics model + new single diffractive
 */
class ScatteringModelSixTrackSD: public ScatteringModelFixed
{
public:
	ScatteringModelSixTrackSD(MaterialProperties* mat);
};

}
#endif
