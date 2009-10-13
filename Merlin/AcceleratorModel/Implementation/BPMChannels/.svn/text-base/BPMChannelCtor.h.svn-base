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

#ifndef BPMChannelCtor_h
#define BPMChannelCtor_h 1

#include "merlin_config.h"
// ChannelServer
#include "AcceleratorModel/Implementation/ChannelServer.h"
// BPMDataBufferServer
#include "AcceleratorModel/Implementation/BPMChannels/BPMDataBufferServer.h"

class BPMChannel;

//	Specialised constructor for BPM channels. The type of
//	channel constructed depends on the value of the
//	inherited key string (='X' or 'Y').

class BPMChannelCtor : public ChannelServer::ChannelCtor
{
public:

    // Constructor taking either 'x' or 'y' to designate
    // which plane of BPM channel to construct.
    BPMChannelCtor (char xy);

    //	Constructs and returns a BPMChannel.
    virtual ROChannel* ConstructRO (ModelElement* anElement);

    //	Returns NULL.
    virtual RWChannel* ConstructRW (ModelElement* anElement);

    static BPMDataBufferServer theServer;
};


#endif
