// class QuadIntegrator
// ----------------------------------------------------------
//
// A quadrupole integrator which calculates the exact
// linear map for each particle energy.
//
// This local integrator is intended to override the internal
// particle integrator for multipoles, which using an second-
// order expansion for the energy.

#include "BasicTransport/TransportMatrix.h"
#include "BasicTransport/MatrixMaps.h"
#include "NumericalUtils/PhysicalConstants.h"
#include "QuadIntegrator.h"

using namespace PhysicalUnits;
using namespace PhysicalConstants;

namespace ParticleTracking {
    
    void QuadIntegrator::TrackStep(double ds)
    {
        double dBdx = currentComponent->GetFieldStrength();
        double p0 = currentBunch->GetReferenceMomentum();
        double brho = p0/eV/SpeedOfLight;
        
        RMtrx Rm(2);

        for(ParticleBunch::iterator p = currentBunch->begin(); p!= currentBunch->end(); p++) {
            double k1 = dBdx/(brho*(1+(*p).dp()));
            TransportMatrix::QuadrupoleR(ds,k1,Rm.R);
            Rm.Apply(*p);
        }
    };
    
}; // end namespace ParticleTracking
