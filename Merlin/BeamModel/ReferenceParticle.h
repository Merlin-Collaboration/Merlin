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

#ifndef ReferenceParticle_h
#define ReferenceParticle_h 1

#include "merlin_config.h"
#include <cassert>

//	A ReferenceParticle represents that particle which sits
//	on the nominal orbit. It is responsible for maintaining
//	the reference momentum and time (ct) for the bunch or
//	map. ReferenceParticle cannot be instantiated, but is
//	designed as a mixin for bunch or map-like classes.

class ReferenceParticle
{
public:

    //	Returns the reference momentum in GeV/c.
    double GetReferenceMomentum () const;

    //	Returns the reference time in ct (meters).
    double GetReferenceTime () const;

    //	Returns either +1, 0 or -1.
    double GetChargeSign () const;

    //	Sets the reference momentum to p GeV/c. p must be
    //	greater than zero.
    void SetReferenceMomentum (double p);

    //	Increments the reference momentum by dp GeV/c, returning
    //	the new value.
    double IncrReferenceMomentum (double dp);

    //	Sets the reference time in ct (meters).
    void SetReferenceTime (double ct);

    //	Increments the reference time by dct meters.
    double IncrReferenceTime (double dct);

protected:

    ReferenceParticle (double p, double q = 1);

    //	Sets the charge sign.
    void SetChargeSign (double q);

    // Data Members for Class Attributes

    //	reference momentum in GeV/c
    double p0;

    //	reference time in ct (meters)
    double ct0;

    //	The charge sign of the particles.
    double qs;
};

inline ReferenceParticle::ReferenceParticle (double p, double q)
        : p0(p),ct0(0)
{
    assert(p>0);
    SetChargeSign(q);
}

inline double ReferenceParticle::GetReferenceMomentum () const
{
    return p0;
}

inline double ReferenceParticle::GetReferenceTime () const
{
    return ct0;
}

inline double ReferenceParticle::GetChargeSign () const
{
    return qs;
}

inline void ReferenceParticle::SetReferenceMomentum (double p)
{
    p0=p;
    assert(p0>0);
}

inline double ReferenceParticle::IncrReferenceMomentum (double dp)
{
    p0+=dp;
    assert(p0>0);
    return p0;
}

inline void ReferenceParticle::SetReferenceTime (double ct)
{
    ct0=ct;
}

inline double ReferenceParticle::IncrReferenceTime (double dct)
{
    ct0+=dct;
    return ct0;
}

inline void ReferenceParticle::SetChargeSign (double q)
{
    if(q<0)
        qs=-1;
    else if(q==0)
        qs=0;
    else
        qs=1;
}

#endif
