#ifndef MonitorProcess_h 
#define MonitorProcess_h 1

#include "BeamDynamics/BunchProcess.h"
#include <vector>
#include <string>
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"


class MonitorProcess : public BunchProcess
{
public:

        vector<string> dump_at_elements;
        Bunch* bunch;
        string file_prefix;

        MonitorProcess(const string& aID = "MONITOR",  int prio=0, const string& prefix = "");
        void SetPrefix(const string& prefix);
        void AddElement(const string e);
        void InitialiseProcess (Bunch& this_bunch);
        void DoProcess (const double ds);
        double GetMaxAllowedStepSize() const;
        void SetCurrentComponent (AcceleratorComponent& component);

};

#endif
