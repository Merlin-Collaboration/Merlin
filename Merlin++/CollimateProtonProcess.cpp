/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iostream>
#include <unistd.h>
#include <vector>

#include "merlin_config.h"

#include "Aperture.h"
#include "CollimatorAperture.h"
#include "InterpolatedApertures.h"
#include "Collimator.h"

#include "PSvector.h"
#include "ParticleComponentTracker.h"
#include "ParticleBunch.h"

#include "CollimateProtonProcess.h"
#include "ScatteringProcess.h"
#include "ScatteringModel.h"

#include "utils.h"
#include "PhysicalUnits.h"
#include "PhysicalConstants.h"

using namespace std;
using namespace Collimation;

namespace ParticleTracking
{

CollimateProtonProcess::CollimateProtonProcess(int priority, int mode, std::ostream* osp) :
	CollimateParticleProcess(priority, mode, osp), scattermodel(nullptr)
{
}

CollimateParticleProcess::ScatterOutcome CollimateProtonProcess::DoScatter(Particle& p)
{
	double P0 = currentBunch->GetReferenceMomentum();
	double E0 = sqrt(P0 * P0 + pow(PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV, 2));

	bool scatter_plot = 0;
	bool jaw_impact = 0;

	// Length of the collimator
	double coll_length = currentComponent->GetLength();

	double z = int_s;
	double lengthtogo = s - z;
	Collimator* C = static_cast<Collimator*>(currentComponent);

	string ColName = currentComponent->GetName();
/*
    if(scattermodel->ScatterPlot_on)
    {
        for(vector<string>::iterator its = scattermodel->ScatterPlotNames.begin(); its !=
            scattermodel->ScatterPlotNames.end(); ++its)
        {
            if(ColName == *its)
            {
                scatter_plot = 1;
            }
        }
    }

    if(scattermodel->JawImpact_on)
    {
        for(vector<string>::iterator its = scattermodel->JawImpactNames.begin(); its !=
            scattermodel->JawImpactNames.end(); ++its)
        {
            if(ColName == *its)
            {
                jaw_impact = 1;
            }
        }
    }
 */
	const Aperture *colap = C->GetAperture();

	//set scattering model
	if(scattermodel == nullptr)
	{
		std::cout << "\nCollimateProtonProcess::SoScatter::WARNING: no ScatteringModel set." << std::endl;
		std::cout << "Use 'myCollimateProcess->SetScatteringModel(myScatter);'" << std::endl;
		exit(EXIT_FAILURE);
	}
	while(lengthtogo > 0)
	{
		double E1 = E0 * (1 + p.dp());
		//Note that pathlength should be calculated with E0

		double xlen = scattermodel->PathLength(C->GetMaterialProperties(), E0);

		double E2 = 0;

		bool interacted = (lengthtogo > xlen);
		double step_size = interacted ? xlen : lengthtogo;

		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());

		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();
/*
        //Jaw Impact
        if(jaw_impact && z == 0)
        {
            scattermodel->JawImpact(p, ColParProTurn, ColName);
        }

        //Scatter Plot
        if(scatter_plot && z == 0)
        {
            scattermodel->ScatterPlot(p, z, ColParProTurn, ColName);
        }
 */
		//Energy Loss
		scattermodel->EnergyLoss(p, step_size, C->GetMaterialProperties(), E0);
		E2 = E0 * (1 + p.dp());

		if(E2 <= 1.0)
		{
			p.ct() = z;

			if(CollimationOutputSet)
			{
				for(CollimationOutputIterator = CollimationOutputVector.begin(); CollimationOutputIterator !=
					CollimationOutputVector.end(); ++CollimationOutputIterator)
				{
					(*CollimationOutputIterator)->Dispose(*currentComponent, (z + zstep), p, ColParProTurn);
				}
			}
			return ScatterOutcome::absorbed;
		}

		//MCS
		scattermodel->Straggle(p, step_size, C->GetMaterialProperties(), E1, E2);

		if((E2 < (E0 / 100.0)))
		{
			return ScatterOutcome::absorbed;
		}

		//Check if (returned to aperture) OR (travelled through length)
		z += zstep;
		if(scatter_plot)
		{
			scattermodel->ScatterPlot(p, z, ColParProTurn, ColName);
		}

		if((colap->CheckWithinApertureBoundaries((p.x()), (p.y()), z)))
		{
			// check it does not come back in
			double extrax = p.x() + p.xp() * lengthtogo;
			double extray = p.y() + p.yp() * lengthtogo;
			if(colap->CheckWithinApertureBoundaries(extrax, extray, z + lengthtogo))
			{
				//escaped jaw, so propagate to end of element
				p.x() = extrax;
				p.y() = extray;
				return ScatterOutcome::survived;
			}
		}

		if(xlen > lengthtogo)
		{
			return ScatterOutcome::survived;
		}

		//Scattering - use E2

		if(interacted)
		{
			if(!scattermodel->ParticleScatter(p, C->GetMaterialProperties(), E2))
			{ // lost
				p.ct() = z;

				if(CollimationOutputSet)
				{
					for(CollimationOutputIterator = CollimationOutputVector.begin(); CollimationOutputIterator !=
						CollimationOutputVector.end(); ++CollimationOutputIterator)
					{
						(*CollimationOutputIterator)->Dispose(*currentComponent, (z + zstep), p, ColParProTurn);
					}
				}
				return ScatterOutcome::absorbed;
			} // else {cout<<" CHECK6  scatter new value "<<p.x()<<" "<<p.y()<<" "<<p.ct()<<" "<<p.dp()<<" "<<p<<endl;}
		}

		if((p.dp() < -0.95) || (p.dp() < -1))
		{
			p.ct() = z;

			if(CollimationOutputSet)
			{
				for(CollimationOutputIterator = CollimationOutputVector.begin(); CollimationOutputIterator !=
					CollimationOutputVector.end(); ++CollimationOutputIterator)
				{
					(*CollimationOutputIterator)->Dispose(*currentComponent, (z + zstep), p, ColParProTurn);
				}
			}
			return ScatterOutcome::absorbed;
		}

		lengthtogo -= step_size;

	}

	//Only reach here in the lengthtogo == zero case, otherwise caught in loop
	return ScatterOutcome::survived;
}

void CollimateProtonProcess::SetScatteringModel(Collimation::ScatteringModel* s)
{
	scattermodel = s;
}

} // end namespace ParticleTracking
