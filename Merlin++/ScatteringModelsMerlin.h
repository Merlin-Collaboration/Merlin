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
 * Merlin physics model
 *
 * ScatteringModel preset with Merlin physics model as described in
 * Appleby, R.B., Barlow, R.J., Molson, J.G. et al.
 * Eur. Phys. J. C (2016) 76: 520.
 * https://dx.doi.org/10.1140/epjc/s10052-016-4363-7
 */
class ScatteringModelMerlin: public ScatteringModel
{
public:
	ScatteringModelMerlin() :
		ScatteringModel()
	{
	}
	void Configure(MaterialProperties *, double Energy) override;
};

/**
 * Sixtrack style physics model
 *
 * ScatteringModel preset with physics based on SixTrack K2 scattering
 */
class ScatteringModelSixTrack: public ScatteringModel
{
public:
	ScatteringModelSixTrack() :
		ScatteringModel()
	{
	}
	void Configure(MaterialProperties *, double Energy) override;
};

/**
 * Sixtrack style physics model + new  Ionisation
 */
class ScatteringModelSixTrackIoniz: public ScatteringModel
{
public:
	ScatteringModelSixTrackIoniz() :
		ScatteringModel()
	{
	}
	void Configure(MaterialProperties *, double Energy) override;
};

/**
 * Sixtrack style physics model + new elastic scattering
 */
class ScatteringModelSixTrackElastic: public ScatteringModel
{
public:
	ScatteringModelSixTrackElastic() :
		ScatteringModel()
	{
	}
	void Configure(MaterialProperties *, double Energy) override;
};

/**
 * Sixtrack style physics model + new single diffractive
 */
class ScatteringModelSixTrackSD: public ScatteringModel
{
public:
	ScatteringModelSixTrackSD() :
		ScatteringModel()
	{
	}
	void Configure(MaterialProperties *, double Energy) override;
};

}
#endif
