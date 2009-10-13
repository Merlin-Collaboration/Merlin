#ifndef _h_CouplerWakeFieldProcess
#define _h_CouplerWakeFieldProcess

#include "BeamDynamics/ParticleTracking/WakeFieldProcess.h"
#include "utility/StringPattern.h"
#include <vector>

#include "AcceleratorModel/CombinedWakeRF.h"

namespace ParticleTracking {
//
// a particle wakefield process modified for couple wake and rf kicks
// DK 30 May 2008 - see EUROTeV-2008-003
//
class CouplerWakeFieldProcess : public WakeFieldProcess {
public:

    CouplerWakeFieldProcess(int prio, size_t nb =100, double ns = 3.0);

    virtual void SetCurrentComponent (AcceleratorComponent& component);
    virtual void CalculateWakeT();

private:
    CombinedWakeRF* currentWake;    
    double phi,V,k;
};

}; // end namespace ParticleTracking

#endif
