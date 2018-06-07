/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <fstream>
#include <iomanip>

#include "ScatteringProcess.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

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

void ScatterStuff(PSvector& p, double t, double E0)
{
	double E1 = (p.dp() + 1) * E0;
	double theta = sqrt(t) / E1;
	double phi = RandomNG::uniform(-pi, pi);
	p.xp() += theta * cos(phi);
	p.yp() += theta * sin(phi);
}

void ScatterStuff(PSvector& p, double t, double m, double E0)  // scatter PSvector by t off target m
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

void ScatterStuff(double dp, PSvector& p, double t, double E0)
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
void Rutherford::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	tmin = 0.9982E-3; // DeMolaize thesis page 29 [GeV^2]
	sigma = cs->Get_sig_R();
	E0 = cs->Get_E0();
}

bool Rutherford::Scatter(PSvector& p, double E)
{
	double TargetMass = AtomicMassUnit * mat->GetAtomicMass();

	t = tmin / (1 - RandomNG::uniform(0, 1));
	ScatterStuff(p, t, TargetMass, E0);
	p.type() = 6;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// ST Rutherford
void SixTrackRutherford::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	tmin = 0.9982E-3; // DeMolaize thesis page 29 [GeV^2]
	sigma = cs->Get_sig_R();
	E0 = cs->Get_E0();
}

bool SixTrackRutherford::Scatter(PSvector& p, double E)
{

	t = tmin / (1 - RandomNG::uniform(0, 1));
	ScatterStuff(p, t, E0);
	p.type() = 6;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// Elastic pn
void Elasticpn::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_el();
	E0 = cs->Get_E0();
}
bool Elasticpn::Scatter(PSvector& p, double E)
{
	t = cs->GetElasticScatter()->SelectT();

	ScatterStuff(p, t, AtomicMassUnit, E0);
	p.type() = 3;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// ST Elasticpn
void SixTrackElasticpn::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_el();
	E0 = cs->Get_E0();
}
bool SixTrackElasticpn::Scatter(PSvector& p, double E)
{
	double com_sqd = 2 * ProtonMassMeV * MeV * E;   //ecmsq in SixTrack
	b_pp = 8.5 + 1.086 * log(sqrt(com_sqd)); // slope given on GeV units
	t = -log(RandomNG::uniform(0, 1)) / b_pp;

	ScatterStuff(p, t, E0);
	p.type() = 3;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// Elastic pN
void ElasticpN::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pN_el();
	double b_N_ref = matin->GetSixtrackNuclearSlope();
	b_N = b_N_ref * (cs->Get_sig_pN_tot() / cs->Get_sig_pN_tot_ref());
	E0 = cs->Get_E0();
}
bool ElasticpN::Scatter(PSvector& p, double E)
{
	double TargetMass = AtomicMassUnit * mat->GetAtomicMass();

	t = -log(RandomNG::uniform(0, 1)) / b_N;
	ScatterStuff(p, t, TargetMass, E0);
	p.type() = 2;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// ST Elastic pN
void SixTrackElasticpN::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pN_el();
	double b_N_ref = matin->GetSixtrackNuclearSlope();
	b_N = b_N_ref * (cs->Get_sig_pN_tot() / cs->Get_sig_pN_tot_ref());
	E0 = cs->Get_E0();
}
bool SixTrackElasticpN::Scatter(PSvector& p, double E)
{

	t = -log(RandomNG::uniform(0, 1)) / b_N;
	ScatterStuff(p, t, E0);
	p.type() = 2;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// SD
void SingleDiffractive::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_sd();
	E0 = cs->Get_E0();
}
bool SingleDiffractive::Scatter(PSvector& p, double E)
{
	std::pair<double, double> TM = cs->GetDiffractiveScatter()->Select();
	t = TM.first;
	m_rec = TM.second;
	double com_sqd = (2 * ProtonMassMeV * MeV * E0) + (2 * ProtonMassMeV * MeV * ProtonMassMeV * MeV);
	double dp = m_rec * m_rec * E / com_sqd;

	ScatterStuff(dp, p, t, E0);
	p.type() = 4;
	p.sd() = 1;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

// ST SD
void SixTrackSingleDiffractive::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pn_sd();
	E0 = cs->Get_E0();
}
bool SixTrackSingleDiffractive::Scatter(PSvector& p, double E)
{
	double com_sqd = 2 * ProtonMassMeV * MeV * E0;  //ecmsq in SixTrack
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
	t = -log(RandomNG::uniform(0, 1)) / b;
	dp = xm2 * E / com_sqd;

	ScatterStuff(dp, p, t, E0);
	p.type() = 4;
	p.sd() = 1;

	double E3 = (1 + p.dp()) * E0;
	if(E3 <= 0.1)
	{
		return false;
	}
	else
	{
		return true;
	}
}

//Inelastic
void Inelastic::Configure(Material* matin, CrossSections* CSin)
{
	ScatteringProcess::Configure(matin, CSin);
	sigma = cs->Get_sig_pN_inel();
	E0 = cs->Get_E0();
}

bool Inelastic::Scatter(PSvector& p, double E)
{
	p.type() = 1;
	return false;
} // Particle is lost
