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

#ifndef BPMChannel_h
#define BPMChannel_h 1

#include "merlin_config.h"
// Channels
#include "Channels/Channels.h"

class BPMDataBuffer;

class BPMChannel : public ROChannel
{
public:

    //	Constructor taking the BPMDataBuffer, and the type of
    //	channel ('X' or 'Y').
    BPMChannel (char XorY, BPMDataBuffer* dataBuff);


    //	Returns the ID of the channel (parameter).
    virtual std::string GetID () const;

    //	Returns the current value of the parameter/attribute
    //	associated with the channel.
    virtual double Read () const;

    // Special function for setting arbitrary offset
    void SetBPMOffset(double v);

private:

    char xy;
    BPMDataBuffer* itsData;
};

#endif
