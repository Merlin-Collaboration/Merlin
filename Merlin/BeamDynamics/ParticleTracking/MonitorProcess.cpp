#include "MonitorProcess.h"
#include "AcceleratorModel/AcceleratorComponent.h"
#include <fstream>
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"

using namespace ParticleTracking;

MonitorProcess::MonitorProcess(const string& aID ,int prio, const string& prefix ): ParticleBunchProcess(aID,prio)
{
	//cout << "MonitorProcess()" << endl;
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
	//cout << "MonitorProcess::InitialiseProcess()" << endl;
}
void MonitorProcess::DoProcess (const double ds)
{
	string filename;
	filename = file_prefix + currentComponent->GetName();

	stringstream f;
	f << file_prefix;
	f <<count;
	filename = f.str();
	count++;
	cout << filename << endl;
#ifndef ENABLE_MPI
	//cout << "MonitorProcess::DoProcess(): " << filename << endl;
	ofstream out_file(filename.c_str());
	currentBunch->Output(out_file);
	out_file.close();
#endif

#ifdef ENABLE_MPI
	currentBunch->gather();
	if(currentBunch->MPI_rank == 0)
	{
		//cout << "MonitorProcess::DoProcess(): " << filename << endl;
		ofstream out_file(filename.c_str());
		currentBunch->Output(out_file);
		out_file.close();
	}
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
	std::vector<string>::iterator result = dump_at_elements.begin();
	//active = false;
	active = true;
	while(result != dump_at_elements.end())
	{
		//cout << (*result) << "\t" << component.GetName() << endl;
		if((*result) == component.GetName())
		{
			active = true;
		}
		/*	if(dump_at_elements.begin() == dump_at_elements.end())
			{
				cout << "one" << endl;
				break;
			}
			*/
		result++;
	}

}
