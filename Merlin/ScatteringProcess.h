/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _h_ScatteringProcess
#define _h_ScatteringProcess 1

#include <iostream>
#include <cmath>

#include "merlin_config.h"

#include "PSvector.h"

#include "Material.h"
#include "DiffractiveScatter.h"
#include "ElasticScatter.h"
#include "CrossSections.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

/*

   Definition of the virtual ScatteringProcess class and
   also several child classes derived from it

   Created RJB 23 October 2012
   Modified HR 07.09.2015

 */

namespace Collimation
{

class ScatteringProcess
{
public:
	double sigma;           /// Integrated cross section for this process

protected:
	double E0;              /// Reference energy
	Material* mat;          /// Material of the collimator being hit
	CrossSections* cs;      /// CrossSections object holding all configured cross sections
	double t;               /// Momentum transfer

public:
	virtual ~ScatteringProcess()
	{
	}
	// The first function must be provided for all child classes, and probably the second as well
	virtual bool Scatter(PSvector& p, double E) = 0;
	virtual void Configure(Material* matin, CrossSections* CSin)
	{
		mat = matin;
		cs = CSin;
	}
	virtual std::string GetProcessType() const
	{
		return "ScatteringProcess";
	}
};

/**
 * Rutherford
 */
class Rutherford: public ScatteringProcess
{
	double tmin;
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "Rutherford";
	}
};

class SixTrackRutherford: public ScatteringProcess
{
	double tmin;
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "SixTrackRutherford";
	}
};

/**
 * Elastic pn
 */
class Elasticpn: public ScatteringProcess
{
	double b_pp; /// slope
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "Elastic_pn";
	}
};

class SixTrackElasticpn: public ScatteringProcess
{
	double b_pp; /// slope
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "SixTrackElasic_pn";
	}
};

/**
 * Elastic pN
 */
class ElasticpN: public ScatteringProcess
{
	double b_N; /// slope
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "Elastic_pN";
	}
};

class SixTrackElasticpN: public ScatteringProcess
{
	double b_N; /// slope
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "SixTrackElasic_pN";
	}
};

/**
 * Single Diffractive
 */
class SingleDiffractive: public ScatteringProcess
{
	double m_rec; /// recoil mass
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "SingleDiffractive";
	}

};

class SixTrackSingleDiffractive: public ScatteringProcess
{
	double m_rec; /// recoil mass
	double dp;
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "SixTrackSingleDiffractive";
	}

};

/**
 * Inelastic
 */
class Inelastic: public ScatteringProcess
{
public:
	void Configure(Material* matin, CrossSections* CSin);
	bool Scatter(PSvector& p, double E);
	std::string GetProcessType() const
	{
		return "Inelastic";
	}
};

} //end namespace Collimation

#endif
