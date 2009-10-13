/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2005/03/29 08:19:43 $
// $Revision: 1.1 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_TransRFIntegrator
#define _h_TransRFIntegrator 1

#include "AcceleratorModel/StdComponent/TransverseRFStructure.h"
#include "BeamDynamics/ParticleTracking/ParticleComponentTracker.h"

namespace ParticleTracking {

class TransRFIntegrator : public ParticleComponentTracker::Integrator<TransverseRFStructure> {
protected:

    void TrackStep(double);
    void TrackEntrance();
    void TrackExit();

private:

    void ApplyEndField(double gsgn);
};

}; // end namespace ParticleTracking

#endif

