/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:53 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_LCAVIntegrator
#define _h_LCAVIntegrator 1

#include "AcceleratorModel/StdComponent/TWRFStructure.h"
#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"

namespace ParticleTracking {

class LCAVIntegrator : public ParticleComponentTracker::Integrator<TWRFStructure> {
protected:

    void TrackStep(double);
    void TrackEntrance();
    void TrackExit();

private:

    void ApplyEndField(double gsgn);
};

}; // end namespace ParticleTracking

#endif

