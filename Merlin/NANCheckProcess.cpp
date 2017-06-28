#include "NANCheckProcess.h"
#include "AcceleratorComponent.h"
#include "ParticleBunchProcess.h"
#include "ParticleBunch.h"
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
	ParticleBunch::iterator p;
	size_t count = 0;
	for(p = currentBunch->begin(); p != currentBunch->end(); p++)
	{
		if(std::isnan(p->x()) || std::isnan(p->y()) ||  std::isnan(p->ct()) || std::isnan(p->xp()) || std::isnan(p->yp()) || std::isnan(p->dp()))
		{
			std::cout << "FIRST NAN entry found in particle " << count << " at " << currentComponent->GetQualifiedName() << "\t";
			std::cout << *p << std::endl;
			active = false;
		}
		count++;
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
