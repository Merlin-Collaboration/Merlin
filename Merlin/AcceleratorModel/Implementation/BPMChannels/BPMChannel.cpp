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
// BPMDataBuffer
#include "AcceleratorModel/Implementation/BPMChannels/BPMDataBuffer.h"
// BPMChannel
#include "AcceleratorModel/Implementation/BPMChannels/BPMChannel.h"

BPMChannel::BPMChannel (char XorY, BPMDataBuffer* dataBuff)
        : xy(XorY),itsData(dataBuff)
{
    assert((xy=='X'||xy=='Y')&&(dataBuff!=0));
}

std::string BPMChannel::GetID () const
{
    return std::string("BPM."+itsData->id+'.'+xy);
}

double BPMChannel::Read () const
{
    if(xy == 'X')
        return itsData->x + itsData->x_off;
    else {
        return itsData->y + itsData->y_off;
    }
}

void BPMChannel::SetBPMOffset(double x)
{
    if(xy == 'X')
        itsData->SetOffsets(x,0);
    else
        itsData->SetOffsets(0,x);
}

