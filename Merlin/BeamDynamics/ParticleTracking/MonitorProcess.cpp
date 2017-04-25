#include "MonitorProcess.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include <fstream>
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"

using namespace ParticleTracking;

MonitorProcess::MonitorProcess(const string& aID ,int prio, const string& prefix ): ParticleBunchProcess(aID,prio)
{
	active = true;
	file_prefix = prefix;
	count = 1;
}

void MonitorProcess::SetPrefix(const string& prefix)
{
	file_prefix = prefix;
}

void MonitorProcess::AddElement(const string e)
{
	dump_at_elements.push_back(e);
}


void MonitorProcess::InitialiseProcess (Bunch& bunch)
{
	ParticleBunchProcess::InitialiseProcess(bunch);
	if(!currentBunch)
	{
		active = false;
	}
}
void MonitorProcess::DoProcess (const double ds)
{
	string filename;
	filename = file_prefix + currentComponent->GetName() + "_" + to_string(count);

	count++;
#ifdef ENABLE_MPI
	currentBunch->gather();
	if(currentBunch->MPI_rank == 0)
#endif
	{
		cout << "MonitorProcess writing" << filename << endl;
		ofstream out_file(filename);
		if(!out_file.good())
		{
			cerr << "Error opening " << filename << endl;
			exit(EXIT_FAILURE);
		}
		currentBunch->Output(out_file);
		out_file.close();
	}
#ifdef ENABLE_MPI
	currentBunch->distribute();
#endif
}

double MonitorProcess::GetMaxAllowedStepSize() const
{
	return 1000;
}

void MonitorProcess::SetCurrentComponent (AcceleratorComponent& component)
{
	currentComponent = &component;
	string name = currentComponent->GetName();
	// active if current component name in dump_at_elements
	active = any_of(dump_at_elements.begin(), dump_at_elements.end(),
	                [&name](string &s){return (s == name);});
}
