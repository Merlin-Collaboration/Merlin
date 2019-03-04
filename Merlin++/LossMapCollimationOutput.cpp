/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "LossMapCollimationOutput.h"

namespace ParticleTracking
{

void LossMapCollimationOutput::Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn)
{
	if(currentComponent != &currcomponent)
	{
		currentComponent = &currcomponent;
	}

	temp.reset();

	temp.ElementName = currentComponent->GetQualifiedName().c_str();
	temp.s = currentComponent->GetComponentLatticePosition();

	// pos is the lost position within the element, position is the exact loss position in the lattice
	temp.position = (pos + temp.s);
	temp.length = currentComponent->GetLength();
	temp.lost = 1;

	//calculate 10cm interval - move to 10cm binning?
	double inter = 0.0;
	bool fin = false;

	do
	{
		if((pos >= inter) && (pos < (inter + 0.1)))
		{
			temp.interval = inter;
			fin = true;
		}
		else
		{
			inter += 0.1;
		}
	} while(fin == false);

	//For Loss Maps
	if(currentComponent->GetType() == "Collimator")
	{
		temp.temperature = LossData::Collimator;
	}
	else
	{
		temp.temperature = LossData::Cold;
	}

	//Now check for warm regions.
	for(std::vector<std::pair<double, double> >::const_iterator WarmRegionsIterator = WarmRegions.begin();
		WarmRegionsIterator != WarmRegions.end(); WarmRegionsIterator++)
	{
		if(currentComponent->GetType() != "Collimator" && temp.s >= WarmRegionsIterator->first && temp.s <=
			WarmRegionsIterator->second)
		{
			temp.temperature = LossData::Warm;
		}
	}

	temp.p = particle;

	//pushback vector
	DeadParticles.push_back(temp);
}

LossMapCollimationOutput::LossMapCollimationOutput(OutputType ot)
{
	otype = ot;
}

void LossMapCollimationOutput::Finalise()
{
	//First sort DeadParticles according to s
	sort(DeadParticles.begin(), DeadParticles.end(), Compare_LossData);

	std::cout << "CollimationOutput:: DeadParticles.size() = " << DeadParticles.size() << std::endl;

	int outit = 0;
	int total = 0;

	switch(otype)
	{
	case nearestelement:
		for(std::vector<LossData>::const_iterator it = DeadParticles.begin(); it != DeadParticles.end(); ++it)
		{
			++total;
			// Start at s = min and push back the first LossData
			if(OutputLosses.size() == 0)
			{
				OutputLosses.push_back(*it);
			}

			// If old element ++loss
			if(it->ElementName == OutputLosses[outit].ElementName)
			{
				OutputLosses[outit].lost += 1;
			}
			// If new element OutputLosses.push_back
			else
			{
				OutputLosses.push_back(*it);
				outit++;
			}
		}
		break;
	case precise:
		for(std::vector<LossData>::const_iterator it = DeadParticles.begin(); it != DeadParticles.end(); ++it)
		{
			++total;
			// Start at s = min and push back the first LossData
			if(OutputLosses.size() == 0)
			{
				OutputLosses.push_back(*it);
			}

			// If position is equal
			if(it->position == OutputLosses[outit].position)
			{
				OutputLosses[outit].lost += 1;
			}
			// If new element outit.push_back
			else
			{
				OutputLosses.push_back(*it);
				outit++;
			}
		}
		break;
	case tencm:
		for(std::vector<LossData>::const_iterator it = DeadParticles.begin(); it != DeadParticles.end(); ++it)
		{
			++total;
			//if no losses are yet stored
			if(OutputLosses.size() == 0)
			{
				OutputLosses.push_back(*it);
				OutputLosses[outit].lost = 1;
			}
			// If in the same bin ++loss
			else if((it->ElementName == OutputLosses[outit].ElementName) && (it->interval ==
				OutputLosses[outit].interval))
			{
				OutputLosses[outit].lost += 1;
			}
			// If new element outit.push_back and set loss to 1
			else
			{
				OutputLosses.push_back(*it);
				outit++;
				OutputLosses[outit].lost = 1;
			}
		}
		break;
	}

	std::cout << "CollimationOutput:: OutputLosses.size() = " << OutputLosses.size() << std::endl;
	std::cout << "CollimationOutput:: Total losses = " << total << std::endl;
}

void LossMapCollimationOutput::Output(std::ostream* os)
{
	switch(otype)
	{
	case nearestelement:
		(*os) << std::setw(34) << std::left << "#Name";
		(*os) << std::setw(34) << std::left << "s";
		(*os) << std::setw(16) << std::left << "loss";
		(*os) << std::setw(16) << std::left << "temperature";
		(*os) << std::setw(16) << std::left << "length" << std::endl;

		for(std::vector<LossData>::const_iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
		{
			(*os) << std::setw(34) << std::left << its->ElementName;
			(*os) << std::setw(34) << std::left << its->s;
			(*os) << std::setw(16) << std::left << its->lost;
			(*os) << std::setw(16) << std::left << its->temperature;
			(*os) << std::setw(16) << std::left << its->length;
			(*os) << std::endl;
		}
		break;

	case precise:
		(*os) << std::setw(34) << std::left << "#Name";
		(*os) << std::setw(34) << std::left << "s";
		(*os) << std::setw(16) << std::left << "position";
		(*os) << std::setw(16) << std::left << "loss";
		(*os) << std::setw(16) << std::left << "temperature" << std::endl;

		for(std::vector<LossData>::const_iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
		{
			(*os) << std::setw(34) << std::left << its->ElementName;
			(*os) << std::setw(34) << std::left << its->s;
			(*os) << std::setw(34) << std::left << its->position;
			(*os) << std::setw(16) << std::left << its->lost;
			(*os) << std::setw(16) << std::left << its->temperature;
			(*os) << std::endl;

		}
		break;

	case tencm:
		(*os) << std::setw(34) << std::left << "#Name";
		(*os) << std::setw(34) << std::left << "s";
		(*os) << std::setw(16) << std::left << "bin_start";
		(*os) << std::setw(16) << std::left << "loss";
		(*os) << std::setw(16) << std::left << "temperature";
		(*os) << std::setw(16) << std::left << "length" << std::endl;

		for(std::vector<LossData>::const_iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
		{
			(*os) << std::setw(34) << std::left << its->ElementName;
			(*os) << std::setw(34) << std::setprecision(15) << std::left << its->s + its->interval;
			(*os) << std::setw(16) << std::setprecision(4) << std::left << its->interval;
			(*os) << std::setw(16) << std::setprecision(15) << std::left << its->lost;
			(*os) << std::setw(16) << std::left << its->temperature;
			(*os) << std::setw(16) << std::left << its->length;
			(*os) << std::endl;

		}
		break;
	}

}

void LossMapCollimationOutput::SetWarmRegion(std::pair<double, double> wr)
{
	WarmRegions.push_back(wr);
}

void LossMapCollimationOutput::ClearWarmRegions()
{
	WarmRegions.clear();
}

std::vector<std::pair<double, double> > LossMapCollimationOutput::GetWarmRegions() const
{
	return WarmRegions;
}

} //End namespace ParticleTracking
