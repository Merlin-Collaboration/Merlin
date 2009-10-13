/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:51 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include <cassert>
// BPM
#include "AcceleratorModel/ActiveMonitors/BPM.h"
// BPMChannelCtor
#include "AcceleratorModel/Implementation/BPMChannels/BPMChannelCtor.h"
// BPMChannel
#include "AcceleratorModel/Implementation/BPMChannels/BPMChannel.h"

BPMDataBufferServer BPMChannelCtor::theServer;

BPMChannelCtor::BPMChannelCtor (char xy)
        : ChannelCtor("BPM",std::string(1,xy))
{
    assert(xy=='X'||xy=='Y');
}

ROChannel* BPMChannelCtor::ConstructRO (ModelElement* anElement)
{
#ifndef NDEBUG
    BPM* bpm = dynamic_cast<BPM*>(anElement);
    assert(bpm!=0);
#else
    BPM* bpm = static_cast<BPM*>(anElement);
#endif

    BPMDataBuffer* dataBuffer = theServer.GetDataBuffer(bpm);
    return new BPMChannel(key[0],dataBuffer);
}

RWChannel* BPMChannelCtor::ConstructRW (ModelElement* anElement)
{
    return (RWChannel*)0;
}

