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

// BPMDataBufferServer
#include "AcceleratorModel/Implementation/BPMChannels/BPMDataBufferServer.h"

BPMDataBuffer* BPMDataBufferServer::GetDataBuffer (BPM* bpm, bool create)
{
    map<BPM*,BPMDataBuffer>::iterator b = dataBuffers.find(bpm);

    if(b!=dataBuffers.end())
        return &(b->second);

    if(!create)
        return (BPMDataBuffer*)0;

    // make new buffer
    BPMDataBuffer* newbuf = &(dataBuffers[bpm]);
    newbuf->SetID(bpm->GetName());
    bpm->AddBuffer(newbuf);
    return newbuf;
}

