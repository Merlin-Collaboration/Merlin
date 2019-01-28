/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "FlukaCollimationOutput.h"

#include "AcceleratorComponent.h"
#include "Collimator.h"
#include "CollimatorAperture.h"

namespace ParticleTracking
{

FlukaCollimationOutput::FlukaCollimationOutput(OutputType ot)
{
	otype = ot;
}

void FlukaCollimationOutput::Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn)
{
	// If current component is a collimator we store the loss, otherwise we do not
	if(currentComponent != &currcomponent)
	{
		currentComponent = &currcomponent;
	}
	temp.reset();

	Collimator* aCollimator = dynamic_cast<Collimator*>(&currcomponent);
	bool is_collimator = aCollimator;

	if(is_collimator)
	{

		temp.ElementName = currentComponent->GetQualifiedName().c_str();
		temp.s = currentComponent->GetComponentLatticePosition();

		// For FLUKA output pos is the lost position within the element
		temp.position = pos;
		temp.length = currentComponent->GetLength();
		temp.lost = 1;

		temp.coll_id = currentComponent->GetCollID();
		temp.turn = turn;

		const CollimatorAperture* tap = dynamic_cast<const CollimatorAperture*> (currentComponent->GetAperture());
		temp.angle = tap->GetCollimatorTilt();

		temp.p = particle;

		//pushback vector
		DeadParticles.push_back(temp);
	}
}

void FlukaCollimationOutput::Finalise()
{
	for(std::vector<LossData>::const_iterator its = DeadParticles.begin(); its != DeadParticles.end(); ++its)
	{
		if(its->p.type() == 1 || its->p.type() == 4)
		{
			OutputLosses.push_back(*its);
		}
	}
}

void FlukaCollimationOutput::Output(std::ostream* os)
{
	std::cout << std::endl << "FlukaCollimationOutput OutputLosses size = " << OutputLosses.size()
			  << ", DeadParticles.size() = " << DeadParticles.size() << std::endl;
	(*os) << "#\t1=icoll\t2=c_rotation\t3=s\t4=x\t5=xp\t6=y\t7=yp\t8=nabs\t9=np\t10=ntu" << std::endl;
	for(std::vector<LossData>::const_iterator its = OutputLosses.begin(); its != OutputLosses.end(); ++its)
	{
		(*os) << std::setw(16) << std::left << its->coll_id;
		(*os) << std::setw(20) << std::left << its->angle;
		(*os) << std::setw(20) << std::left << its->position;
		(*os) << std::setw(20) << std::left << its->p.x();
		(*os) << std::setw(20) << std::left << its->p.xp();
		(*os) << std::setw(20) << std::left << its->p.y();
		(*os) << std::setw(20) << std::left << its->p.yp();
		(*os) << std::setw(20) << std::left << its->p.type();
		(*os) << std::setw(20) << std::left << its->p.id();
		(*os) << std::setw(20) << std::left << its->turn;
		(*os) << std::endl;
	}
}

} //End namespace
