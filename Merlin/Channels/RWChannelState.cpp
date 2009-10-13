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

#include "Channels/RWChannelState.h"

RWChannelState::RWChannelState(const std::vector<RWChannel*>& channels)
        : chs(channels), values(channels.size())
{
    CacheCurrentState();
}

void RWChannelState::Reset()
{
    chs.WriteAll(values);
}

void RWChannelState::CacheCurrentState()
{
    chs.ReadAll(values);
}

