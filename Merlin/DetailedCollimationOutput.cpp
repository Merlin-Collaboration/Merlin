#include "DetailedCollimationOutput.h"

namespace ParticleTracking
{

DetailedCollimationOutput::DetailedCollimationOutput(OutputType ot)
{
	otype = ot;
}

void DetailedCollimationOutput::Dispose(AcceleratorComponent& currcomponent, double pos, Particle& particle, int turn)
{
	if (currentComponent != &currcomponent)
	{
		currentComponent = &currcomponent;
	}
	temp.reset();
	temp.ElementName = currentComponent->GetName();

	bool active = any_of(ids.begin(), ids.end(),
	                     [this](StringPattern &s)
	{
		return (s.Match(this->temp.ElementName));
	});

	if(active)
	{
		temp.s = currentComponent->GetComponentLatticePosition();
		temp.position = pos;
		temp.length = currentComponent->GetLength();
		temp.lost = 1;
		temp.turn = turn;
		temp.p = particle;
		DeadParticles.push_back(temp);
	}
}

void DetailedCollimationOutput::Finalise()
{
}

void DetailedCollimationOutput::Output(std::ostream* os)
{
	(*os) << "#name s pos x xp y yp type id location turn" << std::endl;
	for(auto its = DeadParticles.begin(); its != DeadParticles.end(); ++its)
	{
		(*os) << std::setw(16) << std::left << its->ElementName;
		(*os) << std::setw(20) << std::left << its->s;
		(*os) << std::setw(20) << std::left << its->position;
		(*os) << std::setw(20) << std::left << its->p.x();
		(*os) << std::setw(20) << std::left << its->p.xp();
		(*os) << std::setw(20) << std::left << its->p.y();
		(*os) << std::setw(20) << std::left << its->p.yp();
		(*os) << std::setw(20) << std::left << its->p.type();
		(*os) << std::setw(20) << std::left << its->p.id();
		(*os) << std::setw(20) << std::left << its->p.location();
		(*os) << std::setw(20) << std::left << its->turn;
		(*os) << std::endl;
	}
}

void DetailedCollimationOutput::AddIdentifier(const std::string e)
{
	ids.push_back(e);
}

}//End namespace

