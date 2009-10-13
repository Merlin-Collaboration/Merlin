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

#ifndef RWChannelState_h
#define RWChannelState_h 1

#include "merlin_config.h"
#include <vector>
#include <iostream>
#include "Channels/Channels.h"

class RWChannelState {
public:

    RWChannelState(const std::vector<RWChannel*>& channels);
    void Reset();
    void CacheCurrentState();

private:

    RWChannelArray chs;
    RealVector values;
};

#endif
