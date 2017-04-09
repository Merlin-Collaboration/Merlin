#ifndef MonitorProcess_h
#define MonitorProcess_h 1

#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
#include <vector>
#include <string>
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

namespace ParticleTracking
{
/**
 * A process for recording particle coordinates at elements
 *
 * Can be attached to the tracker to record particle coordinates at all
 * or specific elements to files.
 */
class MonitorProcess : public ParticleBunchProcess
{
private:
	vector<string> dump_at_elements;
	string file_prefix;
	unsigned int count;

public:
	/**
	 * Create a MonitorProcess
	 *
	 * \param aID Process ID
	 * \param prio Process priority
	 * \param prefix Output file name prefix
	 *
	 */
	MonitorProcess(const string& aID = "MONITOR",  int prio=0, const string& prefix = "");

	/// Set the output file name prefix
	void SetPrefix(const string& prefix);

	/// Add element at which to record
	void AddElement(const string e);
	void InitialiseProcess (Bunch&  bunch);
	void DoProcess (const double ds);
	double GetMaxAllowedStepSize() const;
	void SetCurrentComponent (AcceleratorComponent& component);

};

} // end namespace ParticleTracking
#endif
