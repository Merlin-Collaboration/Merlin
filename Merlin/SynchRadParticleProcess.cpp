/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <cmath>
#include <algorithm>
#include "utils.h"

#include "SynchRadParticleProcess.h"
#include "SectorBend.h"
#include "RectMultipole.h"
#include "RandomNG.h"

#include "PhysicalUnits.h"
#include "PhysicalConstants.h"
#include "NumericalConstants.h"

//#define PHOTCONST1 3*PlanckConstant*eV/2/twoPi/ElectronMass
//#define PHOTCONST2 2*ElectronCharge*ElectronCharge*ElectronCharge*FreeSpacePermeability/9/ElectronMass/PlanckConstant

//#define PHOTCONST1 3*PlanckConstant*eV/2/twoPi/ProtonMass
//#define PHOTCONST2 2*ElectronCharge*ElectronCharge*ElectronCharge*FreeSpacePermeability/9/ProtonMass/PlanckConstant

namespace
{

using namespace PhysicalUnits;
using namespace PhysicalConstants;

struct ApplySR
{

	typedef double (*spec_gen)(double);

	const MultipoleField& Bf;
	double dL;
	double P0;
	bool symp;
	spec_gen photgen;

	double meanU;
	size_t n;

	double PHOTCONST1, PHOTCONST2, ParticleMassMeV;

	ApplySR(const MultipoleField& field, double dl, double p0, bool sV, const double PCONST1, const double PCONST2,
		const double MassMeV, spec_gen sg = nullptr) :
		Bf(field), dL(dl), P0(p0), symp(sV), photgen(sg), meanU(0), n(0), PHOTCONST1(PCONST1), PHOTCONST2(PCONST2),
		ParticleMassMeV(MassMeV)
	{
	}

	double MeanEnergyLoss() const
	{
		return meanU / n;
	}

	void operator()(PSvector& v)
	{
		double B = abs(Bf.GetField2D(v.x(), v.y()));
		double g = P0 * (1 + v.dp()) / ParticleMassMeV;
		double uc = PHOTCONST1 * B * g * g;
		double u = 0;

		if(photgen)
		{
			int nphot = static_cast<int>(RandomNG::poisson((PHOTCONST2 * 15. * sqrt(3.) / 8.) * B * dL));
			for(int n = 0; n < nphot; n++)
			{
				u += photgen(uc);
			}
		}
		else
		{
			u = PHOTCONST2 * B * dL * uc;
		}

		meanU += u;

		double& px = v.xp();
		double& py = v.yp();
		double& dp = v.dp();

		if(symp)
		{

			double ks = sqrt((1.0 + dp) * (1.0 + dp) - px * px - py * py);
			double kx = px / ks;
			double ky = py / ks;

			dp -= u / P0;

			double kz = (1.0 + dp) / sqrt(1.0 + kx * kx + ky * ky);

			px = kx * kz;
			py = ky * kz;
		}

		else
		{
			dp -= u / P0;
			px /= (1.0 + u / P0);
			py /= (1.0 + u / P0);
		}
		n++;
	}

};
} // end of anonymous namespace

// Class SynchRadParticleProcess

namespace ParticleTracking
{

SynchRadParticleProcess::PhotonGenerator SynchRadParticleProcess::pgen = HBSpectrumGen;

bool SynchRadParticleProcess::sympVars = false;

SynchRadParticleProcess::SynchRadParticleProcess(int prio, bool q)

	: ParticleBunchProcess("SYNCHROTRON RADIATION", prio), ns(1), incQ(false), adjustEref(true), dsMax(0)

{

	GeneratePhotons(q);

}

void SynchRadParticleProcess::SetCurrentComponent(AcceleratorComponent& component)
{

	PHOTCONST1 = 3 * PlanckConstant * eV / 2 / twoPi / currentBunch->GetParticleMass();
	PHOTCONST2 = 2 * ElectronCharge * ElectronCharge * ElectronCharge * FreeSpacePermeability / 9
		/ currentBunch->GetParticleMass() / PlanckConstant;
	ParticleMassMeV = currentBunch->GetParticleMassMeV() * MeV;

	SectorBend* bend = nullptr;
	RectMultipole* rmult = nullptr;

	if((bend = dynamic_cast<SectorBend*>(&component)))
	{
		currentField = &(bend->GetField());
	}

	else if(incQ && (rmult = dynamic_cast<RectMultipole*>(&component)))
	{
		currentField = &(rmult->GetField());
	}
	else
	{
		currentField = nullptr;
	}

	int ns1 = (ns == 0) ? 1 + component.GetLength() / dsMax : ns;
	dL = component.GetLength() / ns1;
	nk1 = 0;
	intS = 0;

	active = currentField && currentBunch;

}

void SynchRadParticleProcess::DoProcess(double ds)
{

	if(fequal(intS += ds, (nk1 + 1) * dL))
	{
		double E0 = currentBunch->GetReferenceMomentum();
		double meanU = for_each(
			currentBunch->begin(),
			currentBunch->end(),
			ApplySR(*currentField, dL, E0, sympVars, PHOTCONST1, PHOTCONST2, ParticleMassMeV,
			quantum)).MeanEnergyLoss();

		// Finally we adjust the reference of the
		// bunch to reflect the mean energy loss
		// of all the particles
		if(adjustEref)
		{
			currentBunch->AdjustRefMomentum(-meanU / E0);
		}
		nk1++;
	}
	active = nk1 != ns;

}

double SynchRadParticleProcess::GetMaxAllowedStepSize() const
{

	return (nk1 + 1) * dL - intS;

}

void SynchRadParticleProcess::IncludeQuadRadiation(bool quadsr)
{

	incQ = quadsr;

}

void SynchRadParticleProcess::SetNumComponentSteps(int n)
{

	ns = n;
	dsMax = 0;

}

void SynchRadParticleProcess::SetMaxComponentStepSize(double ds_max)
{

	ns = 0;
	dsMax = ds_max;

}

void SynchRadParticleProcess::GeneratePhotons(bool gp)
{

	quantum = gp ? pgen : nullptr;

}
} // end namespace ParticleTracking
