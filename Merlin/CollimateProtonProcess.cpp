/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include <iterator>
#include <iostream>
#include <typeinfo>
#include <unistd.h>
#include <vector>
#include <algorithm>

#include "merlin_config.h"

#include "Aperture.h"
#include "CollimatorAperture.h"
#include "InterpolatedApertures.h"
#include "Collimator.h"

#include "ParticleComponentTracker.h"
#include "ParticleBunch.h"

#include "CollimateProtonProcess.h"
#include "ScatteringProcess.h"
#include "ScatteringModel.h"

#include "utils.h"
#include "PhysicalUnits.h"

using namespace Collimation;

namespace ParticleTracking
{

CollimateProtonProcess::CollimateProtonProcess(int priority, int mode, std::ostream* osp) :
	CollimateParticleProcess(priority, mode, osp), scattermodel(nullptr)
{

}

/**
 * returns true if particle survives, false if it dies
 */
bool CollimateProtonProcess::DoScatter(Particle& p)
{
	double P0 = currentBunch->GetReferenceMomentum();
	double E0 = sqrt(P0 * P0 + pow(PhysicalConstants::ProtonMassMeV * PhysicalUnits::MeV, 2));

	bool scatter_plot = 0;
	bool jaw_impact = 0;

	// Length of the collimator
	double coll_length = currentComponent->GetLength();

	double z = currentBunch->int_s;
	double lengthtogo = s - z;

	Collimator* C = static_cast<Collimator*>(currentComponent);

	string ColName = currentComponent->GetName();

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

		double xlen = scattermodel->PathLength(C->material, E0);

		double E2 = 0;

		bool interacted = (lengthtogo > xlen);
		double step_size = interacted ? xlen : lengthtogo;

		double zstep = step_size * sqrt(1 - p.xp() * p.xp() - p.yp() * p.yp());

		p.x() += step_size * p.xp();
		p.y() += step_size * p.yp();

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

		//Energy Loss
		scattermodel->EnergyLoss(p, step_size, C->material, E0);

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
			return true;
		}

		//MCS
		scattermodel->Straggle(p, step_size, C->material, E1, E2);

		if((E2 < (E0 / 100.0)))
		{
			return false;
		}

		//Check if (returned to aperture) OR (travelled through length)
		z += zstep;
		if(scatter_plot)
		{
			scattermodel->ScatterPlot(p, z, ColParProTurn, ColName);
		}

		if((colap->CheckWithinApertureBoundaries((p.x()), (p.y()), z)))
		{
			//escaped jaw, so propagate to end of element
			p.x() += p.xp() * lengthtogo;
			p.y() += p.yp() * lengthtogo;
			return false;
		}

		if(xlen > lengthtogo)
		{
			return false;
		}

		//Scattering - use E2
		if(interacted)
		{
			if(!scattermodel->ParticleScatter(p, C->material, E2))
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
				return true;
			}
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
			return true;
		}

		lengthtogo -= step_size;

	}

	//If we reached here the particle hits the end of the collimator, and thus survives
	return true;
}

void CollimateProtonProcess::SetScatteringModel(Collimation::ScatteringModel* s)
{
	scattermodel = s;
}

} // end namespace ParticleTracking
