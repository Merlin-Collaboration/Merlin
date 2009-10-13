// class QuadIntegrator
// ----------------------------------------------------------
//
// A quadrupole integrator which calculates the exact
// linear map for each particle energy.
//
// This local integrator is intended to override the internal
// particle integrator for multipoles, which using an second-
// order expansion for the energy.

#ifndef QuadIntegrator_h
#define QuadIntegrator_h

#include "AcceleratorModel/StdComponent/StandardMultipoles.h"
#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"

namespace ParticleTracking {
    
    class QuadIntegrator : public ParticleComponentTracker::Integrator<Quadrupole>
    {
    public:
        void TrackStep(double ds);
    };
    
}; // end namespace ParticleTracking

#endif