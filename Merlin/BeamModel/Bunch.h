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

#ifndef Bunch_h
#define Bunch_h 1

#include "merlin_config.h"
#include <iostream>
// ReferenceParticle
#include "BeamModel/ReferenceParticle.h"
// PSTypes
#include "BeamModel/PSTypes.h"
// Space2D
#include "EuclideanGeometry/Space2D.h"

class Transform3D;
class Histogram;

//	A charged particle bunch. Bunch fields a set of methods
//	for accessing certain macroscopic quantities of a single
//	bunch. An object of class Bunch represents an ensemble
//	of like particle only. Functions for accessing the full
//	six-dimensional phase space, as well as projection on a
//	single phase plane are provided. Concrete classes should
//	implement algorithm specific models of a bunch.

class Bunch : public ReferenceParticle
{
public:

    //	Constructor taking the initial reference momentum in Ge
    //	V/c, and the total charge/e of the bunch.
    Bunch (double p, double q);

    //	virtual destructor.
    virtual ~Bunch ();

    //	Returns the total charge (in units of e).
    virtual double GetTotalCharge () const = 0;

    virtual PSmoments& GetMoments (PSmoments& sigma) const = 0;
    PSmoments GetMoments() const {
        PSmoments S;
        return GetMoments(S);
    }
    virtual PSmoments2D& GetProjectedMoments (PScoord u, PScoord v, PSmoments2D& sigma) const = 0;
    virtual PSvector& GetCentroid (PSvector& p) const = 0;
    virtual Point2D GetProjectedCentroid (PScoord u, PScoord v) const = 0;

    //	Set the reference momentum to the mean (centroid)
    //	momentum of the bunch. Returns the new value in GeV/c.
    virtual double AdjustRefMomentumToMean () = 0;

    //	Set the reference time to the mean (centroid) time of
    //	the bunch. Returns the new value in ct (meters).
    virtual double AdjustRefTimeToMean () = 0;

    //	Output a bunch-model dependent representation to the
    //	specified stream.
    virtual void Output (std::ostream& os) const = 0;

    //	Used to generate a 1-D profile of the bunch projected
    //	onto the specified coordinate. The total area of the
    //	historgram is normalised to unity.
    virtual Histogram& ProjectDistribution (PScoord axis, Histogram& hist) const = 0;

    //	Apply the specified 3D coordinate transformation to the
    //	bunch. Returns true if successful (note that it may not
    //	be possible to apply a general 3D transformation to some
    //	concrete bunch representations).
    virtual bool ApplyTransformation (const Transform3D& t) = 0;
};

inline Bunch::Bunch (double p, double q)
        : ReferenceParticle(p,q)
{}

inline Bunch::~Bunch ()
{
    // nothing to do.
}

#endif
