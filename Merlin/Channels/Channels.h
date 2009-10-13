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

#ifndef Channels_h
#define Channels_h 1

#include "merlin_config.h"
#include <vector>
#include <string>
#include <iostream>

// LinearAlgebra
#include "TLAS/LinearAlgebra.h"

//	A Read-Only Channel representing a constant floating
//	point attribute or parameter.

class ROChannel
{
public:

    virtual ~ROChannel ();

    //	Returns the ID of the channel (parameter).
    virtual std::string GetID () const = 0;

    //	Returns the current value of the parameter/attribute
    //	associated with the channel.
    virtual double Read () const = 0;
};

//	A Read-Write Channel,  representing a single
//	floating-point attribute or parameter that can be
//	accessed (Read) or modified (Write).

class RWChannel : public ROChannel
{
public:

    //	Sets the current value of the parameter/attribute
    //	associated with the channel.
    virtual void Write (double value) = 0;

    //	Increments the current value of the parameter/attribute
    //	associated with the channel. Returns the final
    //	(incremented) value.
    virtual double Increment (double delta);
};

//	A linear array (vector) of ROChannels. On destruction,
//	an ROChannelArray object will destroy its associated
//	channels.

class ROChannelArray
{
public:
    //	Constructor taking a vector of ROChannel objects.
    explicit ROChannelArray (const std::vector<ROChannel*>& chnls);
    ROChannelArray(ROChannelArray& rhs);
    ROChannelArray ();

    ~ROChannelArray ();

    //	Reads the n-th channel.
    double Read (size_t n) const;

    //	Reads all the channels and returns the results in vec.
    void ReadAll (RealVector& vec) const;

    //	Convertion to a RealVector (containing the current
    //	values of the channels).
    operator RealVector () const;

    //	Prints the vector as a two column table. Column 1
    //	contains the ID's of the channels; column 2 contains the
    //	current values.
    void Print (std::ostream& os) const;

    //	Returns the size of the array.
    size_t Size () const;

    // Initialisation
    size_t SetChannels(const std::vector<ROChannel*>& chnls);

    // channel access
    const ROChannel& operator[](size_t n) const {return *(channels[n]);}
    ROChannel& operator[](size_t n) {return *(channels[n]);}

protected:

    //	Protected constructing taking the size of the array.
    explicit ROChannelArray (size_t n);

    void DestroyChannels();

    std::vector<ROChannel*> channels;
};

class RWChannelArray : public ROChannelArray
{
public:

    //	Constructor taking a vector of RWChannels.
    explicit RWChannelArray (const std::vector<RWChannel*>& chnls);
    RWChannelArray (RWChannelArray& rhs);
    RWChannelArray ();

    size_t SetChannels (const std::vector<RWChannel*>& chnls);

    //	Write a value to the n-th channel.
    void Write (size_t n, double value);

    //	Write the same value to all the channels.
    void WriteAll (double value);

    //	Copy the vector to the channels.
    void WriteAll (const RealVector& values);

    //	Increments the value of the n-th channel.
    void Increment (size_t n, double value);

    //	Increments  all the channels by the same value.
    void IncrementAll (double value);

    //	Increments the channels by the corresponding values in
    //	the supplied vector.
    void IncrementAll (const RealVector& values);

    //	Copy assignment from a RealVector.
    void operator = (const RealVector& v);

private:

    RWChannel* RWCh (size_t n);
    const RWChannel* RWCh (size_t n) const;
};

inline ROChannel::~ROChannel ()
{}

inline double RWChannel::Increment (double delta)
{
    double v=Read()+delta;
    Write(v);
    return v;
}

inline ROChannelArray::ROChannelArray (const std::vector<ROChannel*>& chnls)
        : channels(chnls)
{}

inline ROChannelArray::ROChannelArray(ROChannelArray& rhs)
        : channels()
{
    channels.swap(rhs.channels);
}

inline ROChannelArray::ROChannelArray (size_t n)
        : channels(n)
{}

inline ROChannelArray::ROChannelArray ()
{}

inline size_t ROChannelArray::SetChannels(const std::vector<ROChannel*>& chnls)
{
    DestroyChannels();
    channels = chnls;
    return channels.size();
}

inline double ROChannelArray::Read (size_t n) const
{
    return channels[n]->Read();
}

inline size_t ROChannelArray::Size () const
{
    return channels.size();
}

inline void RWChannelArray::Write (size_t n, double value)
{
    RWCh(n)->Write(value);
}

inline void RWChannelArray::Increment (size_t n, double value)
{
    RWCh(n)->Increment(value);
}

inline void RWChannelArray::operator = (const RealVector& v)
{
    WriteAll(v);
}

inline RWChannel* RWChannelArray::RWCh (size_t n)
{
    return static_cast<RWChannel*>(channels[n]);
}

inline const RWChannel* RWChannelArray::RWCh (size_t n) const
{
    return static_cast<const RWChannel*>(channels[n]);
}

#endif
