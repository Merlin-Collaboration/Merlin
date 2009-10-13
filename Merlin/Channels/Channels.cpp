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
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iomanip>
// deleters
#include "stdext/deleters.h"
// Channels
#include "Channels/Channels.h"
#include <algorithm>

ROChannelArray::~ROChannelArray ()
{
    DestroyChannels();
}

void ROChannelArray::DestroyChannels()
{
    for_each(channels.begin(),channels.end(),deleter<ROChannel>());
}

void ROChannelArray::ReadAll (RealVector& vec) const
{
    assert(vec.size()==Size());
    for(size_t i=0; i<Size(); i++)
        vec[i]=channels[i]->Read();
}

ROChannelArray::operator RealVector () const
{
    RealVector v(Size());
    ReadAll(v);
    return v;
}

void ROChannelArray::Print (std::ostream& os) const
{
    vector<string> idlist(Size());
    vector<double> vallist(Size());
    int n=0;
    size_t i;

    for(i=0; i<Size(); i++) {
        idlist[i]=channels[i]->GetID();
        vallist[i]=channels[i]->Read();
        if(n<idlist[i].length())
            n=idlist[i].length();
    }

    for(i=0; i<Size(); i++) {
        os<<std::setw(n)<<left<<idlist[i].c_str()<<" = ";
        os<<vallist[i]<<endl;
    }
}


RWChannelArray::RWChannelArray (const std::vector<RWChannel*>& chnls)
        : ROChannelArray(chnls.size())
{
    std::copy(chnls.begin(),chnls.end(),channels.begin());
}

RWChannelArray::RWChannelArray(RWChannelArray& rhs)
        : ROChannelArray(rhs.Size())
{
    std::copy(rhs.channels.begin(),rhs.channels.end(),channels.begin());
    rhs.channels.clear();
}

RWChannelArray::RWChannelArray ()
{}

size_t RWChannelArray::SetChannels (const std::vector<RWChannel*>& chnls)
{
    DestroyChannels();
    channels.resize(chnls.size());
    std::copy(chnls.begin(),chnls.end(),channels.begin());
    return channels.size();
}

void RWChannelArray::WriteAll (double value)
{
    for(size_t i=0; i<Size(); i++)
        RWCh(i)->Write(value);
}

void RWChannelArray::WriteAll (const RealVector& values)
{
    assert(values.size()==Size());
    for(size_t i=0; i<Size(); i++)
        RWCh(i)->Write(values[i]);
}

void RWChannelArray::IncrementAll (double value)
{
    for(size_t i=0; i<Size(); i++)
        RWCh(i)->Increment(value);
}

void RWChannelArray::IncrementAll (const RealVector& values)
{
    assert(values.size()==Size());
    for(size_t i=0; i<Size(); i++)
        RWCh(i)->Increment(values[i]);
}

