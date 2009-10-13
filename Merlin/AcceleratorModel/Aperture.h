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

#ifndef Aperture_h
#define Aperture_h 1

#include "merlin_config.h"
#include "EuclideanGeometry/Space3D.h"

//	Represents the cross section of the vacuum pipe or other
//	collimating aperture.

class Aperture
{
public:
    virtual ~Aperture ();

    //	Returns true if the point (x,y,z) is within the aperture.
    virtual bool PointInside (double x, double y, double z) const = 0;

    //	Returns true if the point p is within the aperture.
    bool PointInside (const Point3D& p) const;

    //	Returns the radius to the aperture at location z and
    //	angle phi.
    virtual double GetRadiusAt (double phi, double z) const = 0;
};

inline Aperture::~Aperture ()
{}

inline bool Aperture::PointInside (const Point3D& p) const
{
    return PointInside(p.x,p.y,p.z);
}

#endif
