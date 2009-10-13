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

#ifndef BPMDataBufferServer_h
#define BPMDataBufferServer_h 1

#include "merlin_config.h"
#include <map>
// BPMDataBuffer
#include "AcceleratorModel/Implementation/BPMChannels/BPMDataBuffer.h"

//	Singleton instance of a strorage class for BPMChannel
//	Data objects. BPMChannelData are indexed by their
//	associated BPM pointer.

class BPMDataBufferServer
{
public:

    BPMDataBuffer* GetDataBuffer (BPM* bpm, bool create = true);
    void Dump(ostream& os);

private:

    std::map<BPM*, BPMDataBuffer> dataBuffers;
};

#endif
