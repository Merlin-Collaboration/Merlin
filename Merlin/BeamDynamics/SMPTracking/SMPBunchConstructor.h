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

#ifndef SMPBunchConstructor_h
#define SMPBunchConstructor_h 1

#include "merlin_config.h"
// BeamData
#include "BeamModel/BeamData.h"
// BunchConstructor
#include "BeamModel/BunchConstructor.h"
// SMPBunch
#include "BeamDynamics/SMPTracking/SMPBunch.h"

namespace SMPTracking {

class SMPBunchConstructor : public BunchConstructor {
public:

    // Constructor taking the beam definition, the number of z slices
    // (ns) and the number of SMPs per slice. Total number of SMPs
    // generated will be nsm*ns
    SMPBunchConstructor (const BeamData& beam, size_t ns, size_t nsm);

    ~SMPBunchConstructor ();

    virtual Bunch* ConstructBunch (int bunchIndex =0) const;

    //	Returns typed particle bunch.
    SMPBunch* ConstructSMPBunch () const;

private:

    size_t ns,np;
    BeamData beamdat;
    double nSigZ,nSigDP;
};

}; // end namespace SMPTracking

#endif
