
#ifndef RingDeltaTProcess_h
#define RingDeltaTProcess_h 1

#include "merlin_config.h"

// ParticleBunchProcess
#include "BeamDynamics/ParticleTracking/ParticleBunchProcess.h"
// MultipoleField
#include "AcceleratorModel/StdField/MultipoleField.h"

namespace ParticleTracking {

class RingDeltaTProcess : public ParticleBunchProcess
{
public:
    RingDeltaTProcess (int prio);
    virtual void SetCurrentComponent (AcceleratorComponent& component);
    virtual void DoProcess (double ds);
    virtual double GetMaxAllowedStepSize () const;
    void SetBendScale (double bendscale);
protected:
private:
    double scale;
    double dL;
    double intS;
};

}; // end namespace ParticleTracking
#endif
