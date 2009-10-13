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

#ifndef RectangularGeometry_h
#define RectangularGeometry_h 1

#include "merlin_config.h"
#include "AcceleratorModel/StdGeometry/CenteredGeometry.h"

//	Represents a straight line segment. Transformations from
//	points on a RectangularGeometry are pure translations
//	along the z-axis.

class RectangularGeometry : public CenteredGeometry
{
public:

    //	Constructor taking the length of the rectangular
    //	geometry (total z extent).
    RectangularGeometry (double l);

    //	Returns a translation along the z-axis of (s-s0).
    virtual Transform3D GetGeometryTransform (double s0, double s) const throw (BeyondExtent);

    //	Returns a translation along the z-axis of either +l/2 or
    //	-l/2 for the entrance and exit boundary planes
    //	respectively.
    virtual Transform3D GetGeometryTransform (BoundaryPlane p) const;

    //	Returns a translation along the z-axis of +l.
    virtual Transform3D GetTotalGeometryTransform () const;
};

inline RectangularGeometry::RectangularGeometry (double l)
        :CenteredGeometry(l)
{}

#endif
