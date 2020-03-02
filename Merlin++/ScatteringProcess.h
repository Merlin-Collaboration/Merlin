/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ScatteringProcess
#define _h_ScatteringProcess 1

#include "merlin_config.h"
#include "PSvector.h"

#include "MaterialProperties.h"
#include "DiffractiveScatter.h"
#include "ElasticScatter.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

/*

   Definition of the virtual ScatteringProcess class and
   also several child classes derived from it

   Created RJB 23 October 2012
   Modified HR 07.09.2015

   Big clean up RJB 3/2/20
 */

class Material;

namespace Collimation
{
class CrossSections;

class ScatteringProcess
{
public:
	double sigma;           /// Integrated cross section for this process
protected:
	double E0;              /// Reference energy
	MaterialProperties* mat;          /// Material of the collimator being hit

public:
	virtual ~ScatteringProcess()
	{
	}
	virtual bool Scatter(PSvector& p, double E) const = 0;
};

/**
 * Rutherford
 */
class Rutherford: public ScatteringProcess
{
	double tmin = 0.9982E-3; // DeMolaize thesis page 29   [GeV^2];
public:
	bool Scatter(PSvector& p, double E) const;
	Rutherford(MaterialProperties* m)
	{
		mat = m;
	}
};

class SixTrackRutherford: public ScatteringProcess
{
	double tmin = 0.9982E-3; // DeMolaize thesis page 29   [GeV^2];
public:
	bool Scatter(PSvector& p, double E) const;
	SixTrackRutherford()
	{
	}
};

/**
 * Elastic pn
 */
class Elasticpn: public ScatteringProcess
{
public:
	ParticleTracking::ppElasticScatter* calculations;     // point to pomeron etc tables
	Elasticpn(double Energy);     //ctor added RJB
	bool Scatter(PSvector& p, double E) const;
};

class SixTrackElasticpn: public ScatteringProcess
{
public:
	bool Scatter(PSvector& p, double E) const;
};

/**
 * Elastic pN
 */
class ElasticpN: public ScatteringProcess
{
	MaterialProperties* mymat;
public:
	ElasticpN(double Energy, MaterialProperties* m = 0);   // =0 must go
	bool Scatter(PSvector& p, double E) const;
};

class SixTrackElasticpN: public ScatteringProcess
{
	MaterialProperties* mymat;
public:
	SixTrackElasticpN(MaterialProperties* m = 0);   // =0 must go
	bool Scatter(PSvector& p, double E) const;
};

/**
 * Single Diffractive
 */
class SingleDiffractive: public ScatteringProcess
{
public:
	ParticleTracking::ppDiffractiveScatter* calculations;     // point to Regge stuff
	SingleDiffractive(double Energy);     //ctor added RJB
	bool Scatter(PSvector& p, double E) const;

};

class SixTrackSingleDiffractive: public ScatteringProcess
{
public:
	bool Scatter(PSvector& p, double E) const;

};

/**
 * Inelastic
 */
class Inelastic: public ScatteringProcess
{
public:
	bool Scatter(PSvector& p, double E) const;
};

} //end namespace Collimation

#endif
