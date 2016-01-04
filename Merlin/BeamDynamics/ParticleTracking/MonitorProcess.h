#ifndef MonitorProcess_h
#define MonitorProcess_h 1

#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include <vector>
#include <string>
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

namespace ParticleTracking
{
class MonitorProcess : public ParticleBunchProcess
{
public:

	vector<string> dump_at_elements;
	string file_prefix;

	unsigned int count;

	MonitorProcess(const string& aID = "MONITOR",  int prio=0, const string& prefix = "");
	void SetPrefix(const string& prefix);
	void AddElement(const string e);
	void InitialiseProcess (Bunch&  bunch);
	void DoProcess (const double ds);
	double GetMaxAllowedStepSize() const;
	void SetCurrentComponent (AcceleratorComponent& component);

};

} // end namespace ParticleTracking
#endif
