#include "NANCheckProcess.h"
#include "AcceleratorComponent.h"
#include "ParticleBunchProcess.h"
#include "ParticleBunch.h"
#include <string>


using namespace ParticleTracking;

NANCheckProcess::NANCheckProcess(const string& aID,int prio ): ParticleBunchProcess(aID,prio), detailed(0), cull(0), halt(0)
{
	active = true;
}

void NANCheckProcess::InitialiseProcess (Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
}

void NANCheckProcess::DoProcess (const double ds)
{
	ParticleBunch::iterator p;
	size_t count = 0;
	bool do_cull = 0;
	for(auto &p: *currentBunch)
	{
		if (reported.count(p.id()))
		{
			continue;
		}
		if(std::isnan(p.x()) || std::isnan(p.y()) ||  std::isnan(p.ct()) || std::isnan(p.xp()) || std::isnan(p.yp()) || std::isnan(p.dp()))
		{
			std::cout << "NAN entry found in currentBunch["<< count <<"], p.id = " << p.id() << ", at " << currentComponent->GetQualifiedName() << std::endl;
			Report(p.id());
			reported.insert(p.id());
			if(cull)
			{
				do_cull = 1;
			}
			if(halt)
			{
				std::cout << "Halting on NAN coordinate" << std::endl;
				abort();
			}
		}
		count++;
	}
	if(do_cull)
	{
		DoCull();
	}
}

void NANCheckProcess::Report (int id) const
{
	if(detailed)
	{
		auto p_prev = find_if (prev_coords.begin(), prev_coords.end(), [&id](const PSvector p)
		{
			return p.id() == id;
		});
		auto p_start = find_if (start_coords.begin(), start_coords.end(), [&id](const PSvector p)
		{
			return p.id() == id;
		});
		if (p_prev != prev_coords.end())
		{
			std::cout << "prev     " << *p_prev << std::endl;
		}
		std::cout << "start    " << *p_start << std::endl;
	}

	auto p_cur = find_if (currentBunch->begin(), currentBunch->end(), [&id](const PSvector p)
	{
		return p.id() == id;
	});
	std::cout << "current  " << *p_cur << std::endl;
}

void NANCheckProcess::DoCull ()
{
	ParticleBunch* NewBunch=new ParticleBunch(currentBunch->GetReferenceMomentum(),currentBunch->GetTotalCharge()/currentBunch->size());
	NewBunch->reserve(currentBunch->size());
	for(auto &p: *currentBunch)
	{
		if(!(std::isnan(p.x()) || std::isnan(p.y()) ||  std::isnan(p.ct()) || std::isnan(p.xp()) || std::isnan(p.yp()) || std::isnan(p.dp())))
		{
			NewBunch->AddParticle(p);
		}
	}
	currentBunch->clear();
	currentBunch->swap(*NewBunch);
	delete NewBunch;
}

double NANCheckProcess::GetMaxAllowedStepSize() const
{
	return 100000;
}

void NANCheckProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	currentComponent = &component;
	if(detailed)
	{
		prev_coords.clear();
		prev_coords = start_coords;
		start_coords.clear();
		start_coords = currentBunch->GetParticles();
	}
}
