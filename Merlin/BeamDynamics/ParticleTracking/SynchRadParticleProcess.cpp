//   Read the documentation to learn more about C++ code generator
//   versioning.

/*
* Merlin C++ Class Library for Charged Particle Accelerator Simulations
* 
* Class library version 2.0 (1999)
* 
* file Merlin\BeamDynamics\ParticleTracking\SynchRadParticleProcess.cpp
* last modified 09/06/01 01:50:54 PM
*/

/*
* This file is derived from software bearing the following
* restrictions:
*
* MERLIN C++ class library for 
* Charge Particle Accelerator Simulations
* Copyright (c) 2001 by The Merlin Collaboration.
* - ALL RIGHTS RESERVED - 
*
* Permission to use, copy, modify, distribute and sell this
* software and its documentation for any purpose is hereby
* granted without fee, provided that the above copyright notice
* appear in all copies and that both that copyright notice and
* this permission notice appear in supporting documentation.
* No representations about the suitability of this software for
* any purpose is made. It is provided "as is" without express
* or implied warranty.
*/


#include <cmath>
#include <algorithm>
#include "NumericalUtils/utils.h"

// SynchRadParticleProcess
#include "BeamDynamics/ParticleTracking/SynchRadParticleProcess.h"
// SectorBend
#include "AcceleratorModel/StdComponent/SectorBend.h"
// RectMultipole
#include "AcceleratorModel/StdComponent/RectMultipole.h"
// RandomNG
#include "Random/RandomNG.h"

#include "NumericalUtils/PhysicalUnits.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "NumericalUtils/NumericalConstants.h"

#define PHOTCONST1 3*PlanckConstant*eV/2/twoPi/ElectronMass
#define PHOTCONST2 2*ElectronCharge*ElectronCharge*ElectronCharge*FreeSpacePermeability/9/ElectronMass/PlanckConstant

namespace {

using namespace PhysicalUnits;
using namespace PhysicalConstants;

struct ApplySR {

    typedef double (*spec_gen)(double);

    const MultipoleField& Bf;
    double dL;
    double P0;
	bool symp;
    spec_gen photgen;

    double meanU;
    size_t n;

    ApplySR(const MultipoleField& field, double dl, double p0, bool sV, spec_gen sg =0)
            : Bf(field),dL(dl),P0(p0),symp(sV),photgen(sg),meanU(0),n(0) {};

    double MeanEnergyLoss() const { return meanU/n; }

    void operator()(PSvector& v)
    {
        double B  = abs(Bf.GetField2D(v.x(),v.y()));
        double g  = P0 * (1 + v.dp())/(ElectronMassMeV*MeV);
        double uc = PHOTCONST1 * B * g * g;
        double u  = 0;

        if(photgen) {
            int nphot = static_cast<int>(RandomNG::poisson( (PHOTCONST2*15.*sqrt(3.)/8.) * B * dL));
            for(int n=0; n<nphot; n++)
                u += photgen(uc);
        }
        else
            u = PHOTCONST2 * B * dL * uc;

        meanU += u;

		double& px = v.xp();
		double& py = v.yp();
		double& dp = v.dp();

		if(symp) {

			double ks  = sqrt((1.0+dp)*(1.0+dp) - px*px - py*py);
			double kx  = px/ks;
			double ky  = py/ks;

			dp -= u / P0;

			double kz = (1.0 + dp)/sqrt(1.0 + kx*kx + ky*ky);

			px = kx*kz;
			py = ky*kz;

		} else {

			dp -= u / P0;
			px /= (1.0 + u/P0);
			py /= (1.0 + u/P0);

		};

        n++;
    }

};
}; // end of anonymous namespace


// Class SynchRadParticleProcess

namespace ParticleTracking {


SynchRadParticleProcess::PhotonGenerator SynchRadParticleProcess::pgen = HBSpectrumGen;

bool SynchRadParticleProcess::sympVars = false;

SynchRadParticleProcess::SynchRadParticleProcess (int prio, bool q)

        : ParticleBunchProcess("SYNCHROTRON RADIATION",prio),ns(1),incQ(false),adjustEref(true),dsMax(0)

{

    GeneratePhotons(q);

}



void SynchRadParticleProcess::SetCurrentComponent (AcceleratorComponent& component)
{

    SectorBend* bend =0;
    RectMultipole* rmult =0;

    if(bend = dynamic_cast<SectorBend*>(&component))
        currentField =  &(bend->GetField());
    else if(incQ && (rmult = dynamic_cast<RectMultipole*>(&component)))
        currentField = &(rmult->GetField());
    else
        currentField = 0;

    int ns1 = (ns==0) ? 1 + component.GetLength()/dsMax : ns;
    dL = component.GetLength()/ns1;
    nk1=0;
    intS=0;

    active = currentField && currentBunch;

}


void SynchRadParticleProcess::DoProcess (double ds)
{

    if(fequal(intS+=ds,(nk1+1)*dL))
    {
        double E0 = currentBunch->GetReferenceMomentum();
        double meanU = for_each(
                           currentBunch->begin(),
                           currentBunch->end(),
                           ApplySR(*currentField,dL,E0,sympVars,quantum)).MeanEnergyLoss();

        // Finally we adjust the reference of the
        // bunch to reflect the mean energy loss
        // of all the particles
        if(adjustEref)
            currentBunch->AdjustRefMomentum(-meanU/E0);
        nk1++;
    }
    active = nk1!=ns;

}


double SynchRadParticleProcess::GetMaxAllowedStepSize () const
{

    return (nk1+1)*dL-intS;

}


void SynchRadParticleProcess::IncludeQuadRadiation (bool quadsr)
{

    incQ = quadsr;

}


void SynchRadParticleProcess::SetNumComponentSteps (int n)
{

    ns = n;
    dsMax = 0;

}

void SynchRadParticleProcess::SetMaxComponentStepSize (double ds_max)
{

    ns = 0;
    dsMax = ds_max;

}


void SynchRadParticleProcess::GeneratePhotons (bool gp)
{

    quantum = gp ? pgen : 0;

}
}; // end namespace ParticleTracking

