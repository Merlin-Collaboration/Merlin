/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.4 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef StableOrbits_h
#define StableOrbits_h 1

#include <list>
#include "AcceleratorModel/AcceleratorModel.h"
#include "BeamDynamics/ParticleTracking/ParticleBunch.h"

using namespace std;
using namespace ParticleTracking;

class StableOrbits
{
public:
    StableOrbits(AcceleratorModel* aModel);
    void SelectStable(ParticleBunch& aBunch, list<size_t>* index);

    int SetTurns(int turns);
    int SetObservationPoint(int n);

private:
    AcceleratorModel* theModel;
    int nturns;
    int obspnt;
};

#endif
