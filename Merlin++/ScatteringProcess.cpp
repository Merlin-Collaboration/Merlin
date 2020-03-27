/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <iomanip>

#include "ScatteringProcess.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"
#include "MaterialData.h"

#include "RandomNG.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;
using namespace ParticleTracking;
using namespace Collimation;

/*

   File contains a useful general purpose routine, and then
   two routines (Configure and Scatter) for a whole lot of standard processes
   Note that AtomicMassUnit is in MeV

 */

void ScatterInit(PSvector& p, double t, double E0)
{
	double E1 = (p.dp() + 1) * E0;
	double theta = sqrt(t) / E1;
	double phi = RandomNG::uniform(-pi, pi);
	p.xp() += theta * cos(phi);
	p.yp() += theta * sin(phi);
}

void ScatterInit(PSvector& p, double t, double m, double E0)  // scatter PSvector by t off target m
{
// for delta in GeV (m[GeV] t[GeV^2])
	double delta = t / (2 * m);
	double E1 = (p.dp() + 1) * E0;
	double E2 = E1 - delta;
	p.dp() = (E2 - E0) / E0;
	double theta = sqrt(t) / E2;
	double phi = RandomNG::uniform(-pi, pi);
	p.xp() += theta * cos(phi);
	p.yp() += theta * sin(phi);
}

void ScatterInit(double dp, PSvector& p, double t, double E0)
{
// for diffractive process where delta is calculated in the scatter function
	double E1 = (p.dp() + 1) * E0;
	double E2 = E1 - dp;
	p.dp() = (E2 - E0) / E0;
	double theta = sqrt(t) / E2;
	double phi = RandomNG::uniform(-pi, pi);
	p.xp() += theta * cos(phi);
	p.yp() += theta * sin(phi);
}

// Rutherford
Rutherford::Rutherford(MaterialProperties* m)
{
	scatterType = RutherfordScattering;
	mat = m;
}

bool Rutherford::Scatter(PSvector& p, double E) const
{
	double TargetMass = AtomicMassUnit * mat->A_R(); //  mat->GetAtomicMass();
	double t = tmin / (1 - RandomNG::uniform(0, 1));
	ScatterInit(p, t, TargetMass, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string Rutherford::GetScatterTypeString() const
{
	return "Rutherford";
}

// SixTrack Rutherford
SixTrackRutherford::SixTrackRutherford()
{
	scatterType = SixTrackRutherfordScattering;
}

bool SixTrackRutherford::Scatter(PSvector& p, double E) const
{
	double t = tmin / (1 - RandomNG::uniform(0, 1));
	ScatterInit(p, t, E);

	// scatters off proton not nucleus. Quite possibly wrong, but kept for not for consistency
	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string SixTrackRutherford::GetScatterTypeString() const
{
	return "SixTrackRutherford";
}

// Elastic pn
Elasticpn::Elasticpn(double Energy)
{
	scatterType = ElasticpnScattering;

	// Do pomeron physics for elastic scattering
	std::cout << " creating ppElasticScatter\n";
	calculations = new ParticleTracking::ppElasticScatter(); // CHECK NEED DELETE
	calculations->SetTMin(1e-4); // CHECK want to INCORPORATE ALL These
	calculations->SetTMax(1.00);
	calculations->SetStepSize(1e-4);
	calculations->GenerateTDistribution(Energy);
	sigma = calculations->GetElasticCrossSectionN();
	std::cout << " Elastic cross section " << sigma << std::endl;
}

bool Elasticpn::Scatter(PSvector& p, double E) const
{
	double t = calculations->SelectT();
	// double t = cs->GetElasticScatter()->SelectT();
	ScatterInit(p, t, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string Elasticpn::GetScatterTypeString() const
{
	return "Elasticpn";
}

// SixTrack Elasticpn
SixTrackElasticpn::SixTrackElasticpn()
{
	scatterType = ElasticpNScattering;
}

bool SixTrackElasticpn::Scatter(PSvector& p, double E) const
{
	double com_sqd = 2 * ProtonMassMeV * MeV * E;   //ecmsq in SixTrack
	double b_pp = 8.5 + 1.086 * log(sqrt(com_sqd)); // slope given on GeV units
	double t = -log(RandomNG::uniform(0, 1)) / b_pp;

	ScatterInit(p, t, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string SixTrackElasticpn::GetScatterTypeString() const
{
	return "SixTrackElasticpn";
}

// Elastic pN
ElasticpN::ElasticpN(double E, MaterialProperties* mat)
{
	scatterType = ElasticpNScattering;
	mymat = mat;
}

bool ElasticpN::Scatter(PSvector& p, double E) const
{
	double TargetMass =  mymat->A_H(); // GetAtomicMass();
	double b_N = 14.1 * pow(TargetMass, 0.66); // deMolaize Equation 1.31
	double t = -log(RandomNG::uniform(0, 1)) / b_N;
	ScatterInit(p, t, AtomicMassUnit * TargetMass, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string ElasticpN::GetScatterTypeString() const
{
	return "ElasticpN";
}

SixTrackElasticpN::SixTrackElasticpN(MaterialProperties* mat)
{
	scatterType = SixTrackElasticpNScattering;
	mymat = mat;
}

// ST Elastic pN
bool SixTrackElasticpN::Scatter(PSvector& p, double E) const
{

	double TargetMass =  mymat->A_H(); // GetAtomicMass();
	double b_N = 14.1 * pow(TargetMass, 0.66); // deMolaize Equation 1.31
	double t = -log(RandomNG::uniform(0, 1)) / b_N;
	ScatterInit(p, t, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string SixTrackElasticpN::GetScatterTypeString() const
{
	return "SixTrackElasticpN";
}

// Single Diffractive
SingleDiffractive::SingleDiffractive(double Energy)
{
	scatterType = SingleDiffractiveScattering;

	// Do pomeron physics for Diffractive scattering
	calculations = new ParticleTracking::ppDiffractiveScatter();
	calculations->SetTMin(0.0001);
	calculations->SetTMax(4.0);
	calculations->SetTStepSize(1e-4);
	calculations->SetXiMin(pow(ProtonMassGeV + 0.135, 2) / (2 * ProtonMassGeV * Energy + Energy * Energy));
	calculations->SetXiMax(0.12);
	calculations->SetXiStepSize(1e-6);
	calculations->GenerateDistribution(Energy);
	sigma = calculations->GetDiffractiveCrossSection();
	std::cout << "CHECK Diffractive cross section " << sigma << std::endl;
}
bool SingleDiffractive::Scatter(PSvector& p, double E) const
{
	std::pair<double, double> TM = calculations->Select();
	double t = TM.first;
	double m_rec = TM.second;
	double com_sqd = (2 * ProtonMassMeV * MeV * E) + (2 * ProtonMassMeV * MeV * ProtonMassMeV * MeV);
	double dp = m_rec * m_rec * E / com_sqd;
	ScatterInit(dp, p, t, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string SingleDiffractive::GetScatterTypeString() const
{
	return "SingleDiffractive";
}

// SixTrack Single Diffractive
SixTrackSingleDiffractive::SixTrackSingleDiffractive()
{
	scatterType = SixTrackSingleDiffractiveScattering;
}

bool SixTrackSingleDiffractive::Scatter(PSvector& p, double E) const
{
	double com_sqd = 2 * ProtonMassMeV * MeV * E;  //ecmsq in SixTrack
	double b_pp = 8.5 + 1.086 * log(sqrt(com_sqd)); // slope given on GeV units
	double xm2 = exp(RandomNG::uniform(0, 1) * log(0.15 * com_sqd));
	double b = 0.0;
	if(xm2 < 2.0)
	{
		b = 2 * b_pp;
	}
	else if(2.0 <= xm2 && xm2 <= 5.0)
	{
		b = (106.0 - 17.0 * xm2) * b_pp / 26.0;
	}
	else if(xm2 > 5.0)
	{
		b = 7.0 * b_pp / 12.0;
	}
	double t = -log(RandomNG::uniform(0, 1)) / b;
	double dp = xm2 * E / com_sqd;

	ScatterInit(dp, p, t, E);

	double E3 = (1 + p.dp()) * E;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

std::string SixTrackSingleDiffractive::GetScatterTypeString() const
{
	return "SixTrackSingleDiffractive";
}


//Inelastic
Inelastic::Inelastic()
{
	scatterType = InelasticScattering;
}

bool Inelastic::Scatter(PSvector& p, double E) const
{
	return false;
}

std::string Inelastic::GetScatterTypeString() const
{
	return "Inelastic";
}
