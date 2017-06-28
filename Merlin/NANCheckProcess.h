#ifndef NANCheckProcess_h
#define NANCheckProcess_h 1

#include "ParticleBunchProcess.h"
#include "ParticleBunch.h"

namespace ParticleTracking
{
class NANCheckProcess : public ParticleBunchProcess
{
public:

	NANCheckProcess(const string& aID = "NAN Check",  int prio=0);
	void InitialiseProcess (Bunch&  bunch);
	void DoProcess (const double ds);
	double GetMaxAllowedStepSize() const;
	void SetCurrentComponent (AcceleratorComponent& component);

};

} // end namespace ParticleTracking
#endif
