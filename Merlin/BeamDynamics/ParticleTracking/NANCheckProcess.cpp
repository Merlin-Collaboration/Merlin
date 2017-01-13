#include "NANCheckProcess.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"
#include <string>

using namespace ParticleTracking;

NANCheckProcess::NANCheckProcess(const string& aID ,int prio ): ParticleBunchProcess(aID,prio) {}

void NANCheckProcess::InitialiseProcess (Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
	active = true;
}

void NANCheckProcess::DoProcess (const double ds)
{
	std::string name = currentComponent->GetQualifiedName();
	ParticleBunch::iterator p;
	for(p = currentBunch->begin(); p != currentBunch->end(); p++)
	{
		if(std::isnan(p->x()) || std::isnan(p->y()) ||  std::isnan(p->ct()) )
		{
			std::cout << "FIRST NAN entry found in " << name << "\t";
			std::cout << *p << std::endl;
			active = false;
		}
	}

}

double NANCheckProcess::GetMaxAllowedStepSize() const
{
	return 100000;
}

void NANCheckProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	currentComponent = &component;
}
