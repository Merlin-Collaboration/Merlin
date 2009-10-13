/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef CenteredGeometry_h
#define CenteredGeometry_h 1

#include "merlin_config.h"
#include "AcceleratorModel/AcceleratorGeometry.h"

//	Base class for accelerator geometry types whose origin
//	is at the center of the their arc length. Hence their
//	extent extends from -length/2 to +length/2.

class CenteredGeometry : public AcceleratorGeometry
{
public:

    //	Constructor taking the arc length of the geometry.
    explicit CenteredGeometry (double l);

    //	Returns the total arc-length of the geometry.
    virtual double GetGeometryLength () const;

    //	Returns the local extent of this geometry.
    virtual AcceleratorGeometry::Extent GetGeometryExtent () const;

protected:

    //	Used to check if (s1s2) is within the geometry bounds.
    void CheckBounds (double s1, double s2) const throw (BeyondExtent);

    //	Used to check if s1 is within the geometry bounds.
    void CheckBounds (double s) const throw (BeyondExtent);

    double len;
};

inline CenteredGeometry::CenteredGeometry (double l)
        : AcceleratorGeometry(),len(l)
{}

inline double CenteredGeometry::GetGeometryLength () const
{
    return len;
}

inline AcceleratorGeometry::Extent CenteredGeometry::GetGeometryExtent () const
{
    double s=len/2;
    return Extent(-s,s);
}

inline void CenteredGeometry::CheckBounds (double s1, double s2) const throw (BeyondExtent)
{
    double hl=len/2;
    if(fabs(s1)>hl || fabs(s2)>hl)
        throw BeyondExtent();
}

inline void CenteredGeometry::CheckBounds (double s) const throw (BeyondExtent)
{
    if(fabs(s)>len/2)
        throw BeyondExtent();
}

#endif
